function [percTimeFix,fixAngle,D] = analyzeVMRgazeData(Exp,D,minMaxPercDistance,plotTrials)
% analyzeVMRgazeData  Analyze gaze data from VMR experiments
%
% [percTimeFix,fixAngle,D] = analyzeVMRgazeData(Exp,D,minMaxPercDistance)
% takes input structs Exp and D created by processVMRtabletData.m and
% detects and analyzes fixations between a minimum and maximum percentage
% of target distance defined in minMaxPercDistance [min% max%]. See below
% for analysis steps.
% Output variables added to D: correctGazeXY, percNoGaze, analyzedGazeData,
% distGaze_meanFix (all fixations), dirGaze_meanFix (only fixations between
% [min% max%] distance).
%
% If plotTrials is set to TRUE, each trial with the detected fixations is
% plotted for visual inspection.
%
% For each trial:
% 1) Correct gaze data (optional) and determine number of useful samples
% 2) Detect saccades and fixations
% 3) Compute and save mean fixation distances and angles
% 4) Compute percentage of time that fixation occurs in each pearl bin
%    during target preview and reach execution
% 5) Create matrix with fixation angles relative to target(s) during target
%    preview and reach execution
% 6) Plot all fixations (optional, plotTrials==true)

% MIT License
% Copyright (c) 2020 Anouk de Brouwer


if nargin==3
    plotTrials = false;
end

fig1 = scaledFigure(2.5,0.9);
states = {'preview','reach'};
load('sacDetParams_VMR.mat')
fs = Exp.setup.fs;

% targets
targetDistance = Exp.targetDistance;
minMaxDistance = [minMaxPercDistance(1)/100*targetDistance minMaxPercDistance(2)/100*targetDistance];
reachTarget = Exp.reachTarget;
target1Angle = Exp.target1Angle;
target2Angle = Exp.target2Angle;
target1Only = ~isnan(target1Angle) & isnan(target2Angle);
target2Only = isnan(target1Angle) & ~isnan(target2Angle);
targetsBoth = ~isnan(target1Angle) & ~isnan(target2Angle);
nTargets = any(~isnan(target1Angle)) + any(~isnan(target2Angle));
nTrials = length(target1Angle);

% pearls
if isfield(Exp,'pearlNumAngle') && all(~isnan(Exp.pearlNumAngle(:,2)))
    pearlAngle = [Exp.pearlNumAngle(:,2)-180; 180];
    dAngle  = pearlAngle(2)-pearlAngle(1);
    pearlBinAngles = [pearlAngle-0.5*dAngle pearlAngle+0.5*dAngle];
else
    pearlAngle = NaN;
    pearlBinAngles = NaN;
    keyboard
end

% target preview
previewDur_s = Exp.timing.previewDur;
if length(previewDur_s)==1
    previewDur_s = ones(nTrials,1)*previewDur_s;
end

% report and error trials
noReport = ~(strcmp(Exp.trialType,'ReachDM_G') | strcmp(Exp.trialType,'ReachDM_I') | strcmp(Exp.trialType,'ReachDM_K'));
if isfield(D,'feedbackAndMTGood')
    feedbackGood = D.feedbackAndMTGood;
else % remove after adapting calcResultsDecisionMakingAndGaze
    keyboard
    feedbackGood = ismember(D.feedbackMessage,'Good');
    MT = D.tTargetGoLeaveRing(:,4)-D.tTargetGoLeaveRing(:,3);
    feedbackGood(MT>0.4) = false;
    % don't remove error trials, because we can use the report of 'Early' and 'Late' reaches
end

% compute median gaze position at target onset
figure; hold on
plot(D.fixXY_iTarget(:,1),D.fixXY_iTarget(:,2),'.')
correctGazeXY = NaN(size(D.fixXY_iTarget));
for b = 1 : max(Exp.blockNo) % take median of each block
    currBlock = Exp.blockNo==b;
    correctGazeXY(currBlock,1) = nanmedian(D.fixXY_iTarget(currBlock,1));
    correctGazeXY(currBlock,2) = nanmedian(D.fixXY_iTarget(currBlock,2));
end
% compute errors
err_raw = nansum(sqrt(D.fixXY_iTarget(:,1).^2+D.fixXY_iTarget(:,2).^2));
fixXY_iTarget_corr = D.fixXY_iTarget-correctGazeXY;
err_corr = nansum(sqrt(fixXY_iTarget_corr(:,1).^2+fixXY_iTarget_corr(:,2).^2));
% plot raw and corrected gaze positions and check if data need to be corrected
plot(fixXY_iTarget_corr(:,1),fixXY_iTarget_corr(:,2),'.')
axis([-100 100 -100 100]); axis equal; horline(0); vertline(0)
xlabel('x (mm)'); ylabel('y (mm)'); legend('raw','corrected')
title(['Error raw = ' num2str(round(err_raw)) ', error corrected = ' num2str(round(err_corr))])
correctGaze = input('Apply correction to gaze data? Yes(1) or no(0): ');
close
if correctGaze
    D.correctGazeXY = correctGazeXY;
else
    D.correctGazeXY = zeros(size(correctGazeXY));
end

% preallocate variables and struct fields
nNotUpdated = NaN(nTrials,1);
nBlink      = NaN(nTrials,1);
D.percNoGaze = NaN(nTrials,1);
D.analyzedGazeData = false(nTrials,1);
fields_perc = {'atStart','outside','atPearl'};
nCol_perc   = [1 1 length(pearlBinAngles)]; nDepth_perc = [1 1 nTargets];
nCol_fix    = 16;
for st = 1 : length(states)
    for f = 1 : length(fields_perc)
        percTimeFix.(states{st}).(fields_perc{f}) = NaN(nTrials,nCol_perc(f),nDepth_perc(f));
    end
    fixAngle.(states{st}).all = NaN(nTrials,nCol_fix,nTargets);
    fixAngle.(states{st}).all_inclStart = NaN(nTrials,nCol_fix,nTargets);
    fixAngle.(states{st}).iOnset_all = NaN(nTrials,nCol_fix,nTargets);
    fixAngle.(states{st}).iOnset_all_inclStart = NaN(nTrials,nCol_fix,nTargets);
    fixAngle.(states{st}).iOffset_all = NaN(nTrials,nCol_fix,nTargets);
    fixAngle.(states{st}).iOffset_all_inclStart = NaN(nTrials,nCol_fix,nTargets);
end
fixAngle.anchor = NaN(nTrials,1); % does not work yet for 2 targets
fixAngle.iAnchor_relReachOn = NaN(nTrials,1); % does not work yet for 2 targets
fixAngle.iTarget_relReachOn = NaN(nTrials,1); % does not work yet for 2 targets
percTimeFix.pearlBinAngles = pearlBinAngles';

%% Loop over trials

for t = 1 : nTrials
    
    % non-reporting trials with good timing
    if noReport(t) && feedbackGood(t)
        
        %% 1) Correct gaze data (optional) and determine number of useful samples
        
        % gaze data
        xyGaze = D.xyGaze{t} - repmat(D.correctGazeXY(t,:),length(D.xyGaze{t}),1); % correct gaze data
        vGaze = D.vGaze{t};
        dirGaze_raw = xy2compassAngle(xyGaze(:,1),xyGaze(:,2),360);
        dirGaze = [dirGaze_raw-target1Angle(t) dirGaze_raw-target2Angle(t)]; % relative to target 1 and target 2
        distGaze = sqrt(xyGaze(:,1).^2 + xyGaze(:,2).^2);
        notUpdated = D.gazeNotUpdated{t};
        blink = D.blink{t};
        noBlink = ~notUpdated & ~blink;
        
        % relevant trial states
        time_tr = D.time{t}-D.time{t}(1);
        if previewDur_s(t)>0
            preview = D.trialState{t}==6; % ShowScene
        else % no preview, use RT
            preview = D.trialState{t}==7; % FillTarget
        end
        reach = false(size(xyGaze,1),1);
        reach(D.iTargetGoLeaveRingEnd(t,3):D.iTargetGoLeaveRingEnd(t,4)) = true;
        previewReach = [preview reach];
        iPreviewStartEnd = [find(preview,1) find(preview,1,'last')];
        iReachStartEnd = [find(reach,1) find(reach,1,'last')];
        nPreviewReach = iReachStartEnd(2)-iPreviewStartEnd(1)+1;
        
        % number of samples in which gaze data is not updated or there is a blink
        % during relevant trial states
        nNotUpdated(t) = sum(notUpdated(iPreviewStartEnd(1):iReachStartEnd(2)));
        nBlink(t) = sum(blink(iPreviewStartEnd(1):iReachStartEnd(2)));
        D.percNoGaze(t) = (nNotUpdated(t)+nBlink(t))/nPreviewReach*100;
        
        %% 2) Detect saccades and fixations
        
        plotOnOff = false;
        % if gaze data is recorded during at least half of the relevant trial states
        if D.percNoGaze(t) < 50
            D.analyzedGazeData(t) = true;
            
            % detect saccades
            xyvpGaze = [xyGaze vGaze noBlink];
            if plotOnOff; figure(fig1); clf; end
            [onsets,offsets] = saccadeOnsetOffset2(xyvpGaze,[1 length(xyGaze)],fs,sacDetParams);
            if plotOnOff
                xlim([0 length(xyvpGaze)])
                vertline(find(diff(preview)),'k'); vertline(find(diff(reach)),'k');
                keyboard; waitforbuttonpress
            end
            noSaccade = true(length(xyGaze),1);
            if onsets>=1
                for i = 1 : length(onsets)
                    noSaccade(onsets(i):offsets(i)) = false;
                end
            end
            
            % detect fixations of 100 ms and longer
            minSamples = 100/(1000/fs);
            [~,iFixStartEnd] = findIntervals([noBlink noSaccade],minSamples);
            fixation = false(length(xyGaze),1);
            if ~isnan(iFixStartEnd)
                for i = 1 : size(iFixStartEnd,1);
                    int = iFixStartEnd(i,1):iFixStartEnd(i,2);
                    fixation(int) = true;
                end
            end
            
            %% 3) Compute mean fixation distances and angles
            
            % compute fixation angles
            [~,~,~,dirGaze_meanFix,distGaze_meanFix] = calcFixationAngles(...
                dirGaze,distGaze,fixation);
            distGood_fix = distGaze_meanFix>minMaxDistance(t,1) & distGaze_meanFix<minMaxDistance(t,2);
            % save mean fixation distances and angles
            D.distGaze_meanFix{t,1} = distGaze_meanFix;
            D.dirGaze_meanFix{t,1} = NaN(size(dirGaze_meanFix));
            D.dirGaze_meanFix{t,1}(distGood_fix,:) = dirGaze_meanFix(distGood_fix,:);
            % mark start fixations
            dirGaze_meanFix_temp = dirGaze_meanFix;
            dirGaze_meanFix_temp(distGaze_meanFix<minMaxDistance(t,1),[~isnan(target1Angle(t)) ~isnan(target2Angle(t))]) ...
                = 999;
            dirGaze_meanFix_temp(distGaze_meanFix>minMaxDistance(t,2)) = NaN;
            
            %% Plot single trials - for example figure in paper
            
            % fig1 = plotSingleTrial(xyGaze,dirGaze_raw,fixation,distGood_fix,Exp,D,t,fig1);
            
            %% Loop over trial states
            
            for st = 1 : size(previewReach,2)
                
                % current state
                currState = previewReach(:,st);
                
                %% 4) Compute percentage of time that gaze is in each pearl bin
                
                if ~isnan(pearlBinAngles(1))
                    % relative to target 1
                    if target1Only(t)
                        [percTimeInBin,percTime_min,percTime_max] = computePercTimeGazeInBin(dirGaze_meanFix(currState,1),...
                            distGaze_meanFix(currState),minMaxDistance(t,:),pearlBinAngles);
                        percTimeFix.(states{st}).atPearl(t,:,1) = percTimeInBin;
                        % relative to target 2
                    elseif target2Only(t)
                        [percTimeInBin,percTime_min,percTime_max] = computePercTimeGazeInBin(dirGaze_meanFix(currState,2),...
                            distGaze_meanFix(currState),minMaxDistance(t,:),pearlBinAngles);
                        percTimeFix.(states{st}).atPearl(t,:,2) = percTimeInBin;
                    elseif targetsBoth(t)
                        [~,iClosest] = min(abs(dirGaze_meanFix(currState,:)),[],2);
                        iClosest(isnan(dirGaze_meanFix(currState,1))) = NaN;
                        % fixations closest to target 1
                        dirGaze_meanFix_T1 = dirGaze_meanFix(currState,1);
                        dirGaze_meanFix_T1(iClosest==2) = NaN;
                        [percTimeInBin_T1,percTime_min,percTime_max] = computePercTimeGazeInBin(dirGaze_meanFix_T1,...
                            distGaze_meanFix(currState),minMaxDistance(t,:),pearlBinAngles);
                        percTimeFix.(states{st}).atPearl(t,:,1) = percTimeInBin_T1*2; % correct because there are 2 targets
                        % fixations closest to target 2
                        dirGaze_meanFix_T2 = dirGaze_meanFix(currState,2);
                        dirGaze_meanFix_T2(iClosest==1) = NaN;
                        [percTimeInBin_T2,percTime_min,percTime_max] = computePercTimeGazeInBin(dirGaze_meanFix_T2,...
                            distGaze_meanFix(currState),minMaxDistance(t,:),pearlBinAngles);
                        percTimeFix.(states{st}).atPearl(t,:,2) = percTimeInBin_T2*2; % correct because there are 2 targets
                        if abs(100 - (sum(percTimeInBin_T1)+sum(percTimeInBin_T2)+percTime_min+percTime_max)) > 0.1
                            disp('Percentage of time doesnt add up to 100')
                            keyboard
                        end
                    end
                    percTimeFix.(states{st}).atStart(t,1) = percTime_min;
                    percTimeFix.(states{st}).outside(t,1) = percTime_max;
                end
                
                %% 5) Create matrix with fixation angles relative to the target(s)
                
                % get mean fixation angles and distances
                dirGaze_meanFix_currState = D.dirGaze_meanFix{t}(currState,:);
                dirGaze_meanFix_currState_temp = dirGaze_meanFix_temp(currState,:);
                if any(~isnan(dirGaze_meanFix_currState(:)));
                    % save mean fixation angles relative to target 1
                    if ~isnan(target1Angle(t))
                        notNan = ~isnan(dirGaze_meanFix_currState(:,1));
                        mDir = unique(dirGaze_meanFix_currState(notNan,1),'stable');
                        fixAngle.(states{st}).all(t,1:length(mDir),1) = mDir';
                        iFixOn = find(diff(isnan(dirGaze_meanFix_currState(:,1)))==-1)+1;
                        if ~isnan(dirGaze_meanFix_currState(1,1))
                            iFixOn = [-1; iFixOn];
                        end
                        iFixOff = find(diff(isnan(dirGaze_meanFix_currState(:,1)))==1)+1;
                        if ~isnan(dirGaze_meanFix_currState(end,1))
                            iFixOff = [iFixOff; 9999];
                        end
                        fixAngle.(states{st}).iOnset_all(t,1:length(mDir),1) = iFixOn;
                        fixAngle.(states{st}).iOffset_all(t,1:length(mDir),1) = iFixOff;
                        % additional variable including start fixations
                        notNan = ~isnan(dirGaze_meanFix_currState_temp(:,1));
                        ddir = dirGaze_meanFix_currState_temp(notNan,1);
                        mDir_temp = ddir([true; diff(ddir)~=0]);
                        fixAngle.(states{st}).all_inclStart(t,1:length(mDir_temp),1) = mDir_temp;
                        fixAngle.(states{st}).iOnset_all_inclStart(t,mDir_temp~=999,1) = iFixOn;
                        fixAngle.(states{st}).iOffset_all_inclStart(t,mDir_temp~=999,1) = iFixOff;
                    end
                    % save mean fixation angles relative to target 2
                    if ~isnan(target2Angle(t))
                        notNan = ~isnan(dirGaze_meanFix_currState(:,2));
                        mDir = unique(dirGaze_meanFix_currState(notNan,2),'stable');
                        fixAngle.(states{st}).all(t,1:length(mDir),2) = mDir';
                        iFixOn = find(diff(isnan(dirGaze_meanFix_currState(:,2)))==-1)+1;
                        if ~isnan(dirGaze_meanFix_currState(1,2))
                            iFixOn = [0; iFixOn];
                        end
                        iFixOff = find(diff(isnan(dirGaze_meanFix_currState(:,2)))==1)+1;
                        if ~isnan(dirGaze_meanFix_currState(end,2))
                            iFixOff = [iFixOff; 9999];
                        end
                        fixAngle.(states{st}).iOnset_all(t,1:length(mDir),2) = iFixOn;
                        fixAngle.(states{st}).iOffset_all(t,1:length(mDir),2) = iFixOff;
                        % additional variable including start fixations
                        notNan = ~isnan(dirGaze_meanFix_currState_temp(:,2));
                        ddir = dirGaze_meanFix_currState_temp(notNan,2);
                        mDir_temp = ddir([true; diff(ddir)~=0]);
                        fixAngle.(states{st}).all_inclStart(t,1:length(mDir_temp),2) = mDir_temp;
                        fixAngle.(states{st}).iOnset_all_inclStart(t,mDir_temp~=999,2) = iFixOn;
                        fixAngle.(states{st}).iOffset_all_inclStart(t,mDir_temp~=999,2) = iFixOff;
                    end
                end
                
            end % end of loop over states
            
            %% 6) Plot all fixations
            
            if plotTrials && (D.percNoGaze(t)<50)
                %%% code needs to be updated
                % figure(fig1);clf
                % iState = [iPreviewStartEnd iReachStartEnd];
                % stateColors = {'k','k','k','k'};
                % fixColors = [repmat({'g--'},1,size(iFixStartEnd,1)) repmat({'r--'},1,size(iFixStartEnd,1))];
                % subplot3([xyGaze(:,1) xyGaze(:,2) dirGaze(:,reachTarget(t))],[0 length(xyGaze)+1],...
                %     [-minMaxDistance(t,2) minMaxDistance(t,2); -minMaxDistance(t,2) minMaxDistance(t,2); -90 90],...
                %     [iState(:)' iFixStartEnd(:)'],[stateColors fixColors],t)
                % horline(0)
                % keyboard %waitforbuttonpress
            end
            
        end
    end % if noReport(t) && feedbackGood(t)
    
end % end of loop over trials

% display number of correct trials
str = sprintf('%d non-reporting trials with good timing and sufficient gaze data',sum(D.percNoGaze<50));
disp(str)

close(fig1);

end


function [percTimeInBin,percTime_min,percTime_max] = computePercTimeGazeInBin(gazeAngle,gazeDistance,minMaxDistance,binAngles)

% make sure gaze angle is between -180 and 180 deg
gazeAngle(gazeAngle<-180) = gazeAngle(gazeAngle<-180)+360;
gazeAngle(gazeAngle>180) = gazeAngle(gazeAngle>180)-360;

% define when gaze distance is good
yesFix = sum(~isnan(gazeDistance));
distanceGood = gazeDistance>minMaxDistance(1) & gazeDistance<minMaxDistance(2);

% count number of samples in each bin and compute percentage of time that
% gaze is in each bin
nBins = size(binAngles,1);
nInBin = NaN(1,nBins);
for b = 1 : nBins
    inBin = gazeAngle>binAngles(b,1) & gazeAngle<binAngles(b,2);
    nInBin(b) = sum(inBin&distanceGood);
end
nYesFix = sum(yesFix);
percTimeInBin = nInBin / nYesFix * 100;

% count number of samples
percTime_min = sum(gazeDistance<minMaxDistance(1)) / nYesFix * 100;
percTime_max = sum(gazeDistance>minMaxDistance(2)) / nYesFix * 100;

end


function [iStartEnd,mAngle,mDist,gazeAngle_meanFix,gazeDistance_meanFix] = calcFixationAngles(gazeAngle,gazeDistance,fixation)

mAngle = [];
mDist = [];
gazeAngle_meanFix = NaN(size(gazeAngle));
gazeDistance_meanFix = NaN(size(gazeAngle,1),1);
[~,iStartEnd] = findIntervals(fixation,10);
% compute mean angle and distance of each fixation
if ~isnan(iStartEnd(1))
    mAngle = NaN(size(iStartEnd,1),size(gazeAngle,2));
    mDist = NaN(size(iStartEnd,1),1);
    for i = 1 : size(iStartEnd,1)
        currFix = iStartEnd(i,1):iStartEnd(i,2);
        % fixation distance
        md = mean(gazeDistance(currFix));
        mDist(i) = md;
        gazeDistance_meanFix(currFix,1) = md;
        % fixation angle
        ga = gazeAngle(currFix,:);
        gar = deg2rad(ga);
        gaur = unwrap(gar); % unwrap gaze angle
        gau = rad2deg(gaur);
        dga = abs(ga-gau);
        mingau = min(gau); maxgau = max(gau);
        mAngle(i,:) = mean(gau);
        if (any((abs(maxgau-mingau))>30) && md>50)
            % check averaging of gaze angle
            plot(iStartEnd(i,1):iStartEnd(i,2),gazeAngle(currFix,1),'b:'); hold on;
            plot(iStartEnd(i,1):iStartEnd(i,2),gau(:,1),'b');
            plot(iStartEnd(i,1):iStartEnd(i,2),gazeAngle(currFix,2),'r:'); hold on;
            plot(iStartEnd(i,1):iStartEnd(i,2),gau(:,2),'r'); hold off
            horline(mAngle(i,1),'b--'); horline(mAngle(i,2),'r--');
            text(iStartEnd(i,1),mAngle(i,1),num2str(mAngle(i,1)),'VerticalAlignment','bottom');
            text(iStartEnd(i,1),mAngle(i,2),num2str(mAngle(i,2)),'VerticalAlignment','bottom');
            xlabel('Time (samples)'); ylabel('Angle (deg)')
            title('Check averaging of gaze angle')
            keyboard;
        end
        mAngle(mAngle<-180) = mAngle(mAngle<-180)+360; % make sure angles are within [-180 180]
        mAngle(mAngle>180) = mAngle(mAngle>180)-360;
        gazeAngle_meanFix(currFix,1) = mAngle(i,1);
        gazeAngle_meanFix(currFix,2) = mAngle(i,2);
        
    end
end

end

function [fig1] = plotSingleTrial(xyGaze,dirGaze_raw,fixation,distGood_fix,Exp,D,t,fig1)

fadedGreen1 = [0 0.5 0]+(1-[0 0.5 0])*0.6;
fadedGreen2 = [0 0.5 0]+(1-[0 0.5 0])*0.8;

% intervals for plotting
time_t = D.time{t}-D.time{t}(1);
int = D.iTargetGoLeaveRingEnd(t,1):D.iTargetGoLeaveRingEnd(t,4);
int_prev = D.iTargetGoLeaveRingEnd(t,1):D.iTargetGoLeaveRingEnd(t,2);
int_RT = D.iTargetGoLeaveRingEnd(t,2):D.iTargetGoLeaveRingEnd(t,3);
int_rch = D.iTargetGoLeaveRingEnd(t,3):D.iTargetGoLeaveRingEnd(t,4);

% pearl angle
pearlAngle = [Exp.pearlNumAngle(:,2)-180; 180];

% hit angle
hitAngle = D.hitAngle_cursor(t,:) + [Exp.target1Angle(t) Exp.target2Angle(t)];
hitAngle = hitAngle(~isnan(hitAngle));
hitXY = Exp.targetDistance(t)*[sind(hitAngle) cosd(hitAngle)];

% counter angles
if Exp.cursorRotation(t)~=0
    counter1Angle = Exp.target1Angle(t) - Exp.target1Rotation(t);
    counter2Angle = Exp.target2Angle(t) - Exp.target2Rotation(t);
    counter1XY = Exp.targetDistance(t)*[sind(counter1Angle) cosd(counter1Angle)];
    counter2XY = Exp.targetDistance(t)*[sind(counter2Angle) cosd(counter2Angle)];
else
    counter1Angle = NaN; counter2Angle = NaN;
    counter1XY = [NaN NaN]; counter2XY = [NaN NaN];
end

figure(fig1); clf;
% plot data in xy coordinates
subplot(1,3,1); hold on
plot(100*sind(pearlAngle),100*cosd(pearlAngle),'ko'); % plot ring
plot(xyGaze(int_prev,1),xyGaze(int_prev,2),'.-','color',[0 0.5 0],'lineWidth',2); % plot gaze
plot(xyGaze(int_RT,1),xyGaze(int_RT,2),'.-','color',fadedGreen1,'lineWidth',2);
plot(xyGaze(int_rch,1),xyGaze(int_rch,2),'.-','color',fadedGreen2,'lineWidth',2);
%plot(D.xyPen{t}(int,1),D.xyPen{t}(int,2),'c-'); % plot hand
plot(D.xyCursor{t}(int,1),D.xyCursor{t}(int,2),'co'); % plot cursor
plot(Exp.target1XY(t,1),Exp.target1XY(t,2),'ro','markerfacecolor','r','markersize',8); % plot target 1
plot([0 Exp.target1XY(t,1)],[0 Exp.target1XY(t,2)],'r');
plot(Exp.target2XY(t,1),Exp.target2XY(t,2),'bo','markerfacecolor','b','markersize',8); % plot target 2
plot([0 Exp.target2XY(t,1)],[0 Exp.target2XY(t,2)],'b')
plot(hitXY(1),hitXY(2),'co','markerfacecolor','c','markersize',8) % plot hit
plot([0 counter1XY(1)],[0 counter1XY(2)],'r--'); % plot aimpoint 1
plot([0 counter2XY(1)],[0 counter2XY(2)],'b--'); % plot aimpoint 2
axis([-125 125 -125 125]); axis equal
xlabel('X (mm)'); ylabel('Y (mm)')
title(['Block ' num2str(Exp.blockNo(t)) ' trial ' num2str(Exp.trialNo(t))])
hold off

% plot xy over time
subplot(1,3,2);
axis([time_t(int(1)) time_t(int(end)) -130 125])
horline(0,'k-')
horline(Exp.target1XY(t,1),'r-'); horline(Exp.target1XY(t,2),'r-.');
horline(Exp.target2XY(t,1),'b-'); horline(Exp.target2XY(t,2),'b-.');
horline(counter1XY,'r--')
horline(counter2XY,'b--')
vertline(D.tTargetGoLeaveRingEnd(t,2:4));
hold on
plot(time_t(fixation),-120,'.','color',[0.5 0.5 0.5])
pg = plot(time_t(int_prev),xyGaze(int_prev,1),'-',... % plot gaze
    time_t(int_prev),xyGaze(int_prev,2),'-.','color',[0 0.5 0],'lineWidth',2);
plot(time_t(int_RT),xyGaze(int_RT,1),'-',...
    time_t(int_RT),xyGaze(int_RT,2),'-.','color',fadedGreen1,'lineWidth',2);
plot(time_t(int_rch),xyGaze(int_rch,1),'-',...
    time_t(int_rch),xyGaze(int_rch,2),'-.','color',fadedGreen2,'lineWidth',2);
plot(time_t([int_RT int_rch]),D.xyCursor{t}([int_RT int_rch],1),'.-',...
    time_t([int_RT int_rch]),D.xyCursor{t}([int_RT int_rch],2),'-.','color','c','linewidth',2); % plot cursor
legend(pg,{'x','y'})
xlabel('Time (s)'); ylabel('Position (mm)')
title('XY positions over time')

% plot angles over time
subplot(1,3,3);
axis([time_t(int(1)) time_t(int(end)) 0 360])
set(gca,'Ytick',0:45:360)
horline(Exp.target1Angle(t),'r-');
horline(Exp.target2Angle(t),'b-');
horline(counter1Angle,'r--')
horline(counter2Angle,'b--')
vertline(D.tTargetGoLeaveRingEnd(t,2:4));
hold on
dirGaze_raw_fixOnly = NaN(size(dirGaze_raw,1),1);
dirGaze_raw_fixOnly(distGood_fix) = dirGaze_raw(distGood_fix);
plot(time_t(int_prev),dirGaze_raw_fixOnly(int_prev),'-','color',[0 0.5 0],'linewidth',2);
plot(time_t(int_RT),dirGaze_raw_fixOnly(int_RT),'-','color',fadedGreen1,'linewidth',2);
plot(time_t(int_rch),dirGaze_raw_fixOnly(int_rch),'-','color',fadedGreen2,'linewidth',2);
dirCursor = xy2compassAngle(D.xyCursor{t}(:,1),D.xyCursor{t}(:,2),360);
plot(time_t(int_rch),dirCursor(int_rch),'-','color','c','linewidth',2);
xlabel('Time (s)'); ylabel('Angle (deg)')
title('Fixation directions over time')

keyboard

end