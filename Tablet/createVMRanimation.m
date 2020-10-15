function createVMRanimation(dataPath,trials,plotGaze)
% createVMRanimation  Create avi animation of VMR task

% MIT License
% Copyright (c) 2020 Anouk de Brouwer

% define path where animation is saved
iSep = strfind(dataPath,filesep);
projectPath = dataPath(1:iSep(end-1)); % two levels up

%% Load data

load(dataPath)
goodTrials = strcmp(D.feedbackMessage(trials),'Good');
if any(~goodTrials)
    disp('Trial selection contains invalid trials.')
    keyboard
end
t = trials(1);

%% Get stimuli
% to do: load raw data file and get sizes and colors from header

fs = 50;
monitorSize = Exp.setup.monitorSize_mm;
startXY = [0 0];
startRadius = Exp.stim.startRadius;
targetDistance = Exp.targetDistance(t);
targetRadius = Exp.stim.targetRadius;
cursorRadius = Exp.stim.cursorRadius;
gazeRadius = cursorRadius/2;
pearlAngle = Exp.pearlNumAngle(:,2);
pearlRadius = 3;
pearlXY = [targetDistance*cosd(pearlAngle) targetDistance*sind(pearlAngle)];
nPearls = length(pearlAngle);

%% Open video file

v = VideoWriter([projectPath 'VMRtrials_' Exp.subjName '.avi']);
v.FrameRate = fs;
open(v);
clear M
ni=0;

%% Create figure

figure; hold on
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2) monitorSize*3]);
set(gca,'position',[0 0 1 1],'units','normalized')
rectangle('Position',[-monitorSize(1)/2 -monitorSize(2)/2 monitorSize],'FaceColor','k'); % background
axis equal; axis([-monitorSize(1)/2 monitorSize(1)/2 -monitorSize(2)/2 monitorSize(2)/2])
set(gca,'Xtick',[]); set(gca,'Ytick',[])

%% Get trial timing and cursor movement

for t = trials
    
    % get target
    targetXY = Exp.target1XY(t,:);
    % get time and cursor movement
    time_tr_ = D.time{t};
    trialState_ = D.trialState{t};
    xyCursor_ = D.xyCursor{t};
    xyGaze_ = D.xyGaze{t};
    hitAngle = -(D.hitAngle_cursor(t)+Exp.target1Angle(t)) + 90; % to unit circle
    xyHit = [cosd(hitAngle) sind(hitAngle)]*targetDistance;
    % define when cursor changed position
    [xy_new,i_new] = unique([xyCursor_ xyGaze_],'rows','stable');
    time_new = time_tr_(i_new)-time_tr_(i_new(1));
    trialState_new = trialState_(i_new);
    % resample
    time_r = (time_new(1):1/fs:time_new(end))'; % multiply fs by two to save at 50% speed
    xyCursor_r = interp1(time_new,xy_new(:,1:2),time_r);
    dCursor_r = sqrt(xyCursor_r(:,1).^2 + xyCursor_r(:,2).^2);
    xyGaze_r = interp1(time_new,xy_new(:,3:4),time_r);
    trialState_r = interp1(time_new,trialState_new,time_r,'nearest');
    nSamples = length(xyCursor_r);
    % get timing of events
    iHit = find(dCursor_r>=targetDistance,1);
    iHit = [iHit-1 iHit];
    dHit = abs(dCursor_r(iHit)-targetDistance);
    iHit = iHit(dHit==min(dHit));
    iTarget = find(trialState_r==6,1);
    iGO = find(trialState_r==7,1);
    if plotGaze && any(any(isnan(xyGaze_r(iTarget:iHit,:))))
        disp('Trial contains blink, discard')
        keyboard
    end
    %% Create animation
    
    % plot pearls, start, and cursor
    for p = 1 : nPearls
        prls(p) = rectangle('Position',[pearlXY(p,1)-pearlRadius...
            pearlXY(p,2)-pearlRadius pearlRadius*2 pearlRadius*2],...
            'Curvature',[1 1],'EdgeColor',[0.8 0.8 0.8]);
    end
    strt = rectangle('Position',[0-startRadius 0-startRadius startRadius*2 startRadius*2],... % start
        'Curvature',[1 1],'EdgeColor','none','FaceColor','w');
    trgt = rectangle('Position',[targetXY-targetRadius targetRadius*2 targetRadius*2],...
        'Curvature',[1 1],'LineWidth',2,'EdgeColor','none','FaceColor','none'); % target
    crsr = rectangle('Position',[0-cursorRadius 0-cursorRadius cursorRadius*2 cursorRadius*2],... % cursor
        'Curvature',[1 1],'EdgeColor',[0 1 1],'FaceColor',[0 1 1]);
    if plotGaze
        gz = rectangle('Position',[0-gazeRadius 0-gazeRadius gazeRadius*2 gazeRadius*2],... % gaze
            'Curvature',[1 1],'EdgeColor','none','FaceColor','y');
    end
    % start trial: display delay
    i = 1;
    while i<iTarget
        crsr.Position = [xyCursor_r(i,:)-cursorRadius cursorRadius*2 cursorRadius*2];
        if plotGaze
            gz.Position = [xyGaze_r(i,:)-gazeRadius gazeRadius*2 gazeRadius*2];
        end
        M(i+ni) = getframe(gcf);
        writeVideo(v,M(i+ni))
        i = i+1;
    end
    % go delay
    trgt.EdgeColor = 'r'; trgt.FaceColor = 'k'; % plot target outline
    while i<iGO
        crsr.Position = [xyCursor_r(i,:)-cursorRadius cursorRadius*2 cursorRadius*2];
        if plotGaze
            gz.Position = [xyGaze_r(i,:)-gazeRadius gazeRadius*2 gazeRadius*2];
        end
        M(i+ni) = getframe(gcf);
        writeVideo(v,M(i+ni))
        i = i+1;
    end
    % reach
    trgt.FaceColor = 'r'; % fill target
    while i<iHit
        crsr.Position = [xyCursor_r(i,:)-cursorRadius cursorRadius*2 cursorRadius*2];
        if plotGaze
            gz.Position = [xyGaze_r(i,:)-gazeRadius gazeRadius*2 gazeRadius*2];
        end
        M(i+ni) = getframe(gcf);
        writeVideo(v,M(i+ni))
        i = i+1;
    end
    % reach complete, freeze cursor for 1 s
    disp(D.hitAngle_cursor(t))
    if abs(D.hitAngle_cursor(t))<4
        trgt.FaceColor = 'g'; trgt.EdgeColor = 'g';
    end
    crsr.FaceColor = 'none'; crsr.LineWidth = 2;
    crsr.Position = [xyHit-cursorRadius cursorRadius*2 cursorRadius*2];
    while i<nSamples
        if plotGaze
            gz.FaceColor = 'y';
            if ~(any(isnan(xyGaze_r(i,:))))
                gz.Position = [xyGaze_r(i,:)-gazeRadius gazeRadius*2 gazeRadius*2];
            end
        end
        M(i+ni) = getframe(gcf);
        writeVideo(v,M(i+ni))
        i = i+1;
    end
    % iti
    prls.delete
    trgt.delete
    strt.delete
    crsr.delete
    if plotGaze
        gz.delete
    end
    while i<(nSamples+fs)
        M(i+ni) = getframe(gcf);
        writeVideo(v,M(i+ni))
        i = i+1;
    end
    M(i+ni) = getframe(gcf);
    writeVideo(v,M(i+ni))
    ni = ni+i;
end

close(v)
close(gcf)
