function processVMRtabletData(projectPath,plotTrials)
% processVMRtabletData  Process raw data from VMR experiments on tablet setup.
%
% processVMRtabletData loads, processes, and visualizes hand and eye
% movement data collected in a visuomotor rotation (VMR) experiment on a
% tablet setup at Queen's University:
% - Wacom tablet with pen and Eyelink 1000 desktop tracker in Abramsky Hall
% - Wacom tablet in mock MRI scanner in Craine Building
% - custom MRI-compatible tablet and Eyelink 1000 remote tracker in fMRI scanner
%
% processVMRtabletData(projectPath) asks the user to specify a project
% folder (if undefined) and select one or multiple participants or sessions.
% The raw data should be stored in a subfolder named '1_RawData', with a
% separate folder containing the raw data (one file per trial) of each
% participant and session. The function loads and cleans the data of each
% trial and saves a single file per participant and session, containing a
% struct 'Exp' with information about the experiment and trials, and a
% struct 'D' with data. The processed data will be saved in a subfolder
% named '2_ProcessedData'.
%
% processVMRtabletData(projectPath,TRUE) plots each trial separately for
% visual inspection during processing, and waits for the user to click a
% button to continue to the next trial. The plot contains a subplot of the
% xy position of the stimuli, hand and gaze, and a subplot of the hand and
% gaze position over time.

% MIT License
% Copyright (c) 2020 Anouk de Brouwer

close all;

if nargin==0
    projectPath = [];
    plotTrials = false;
end

% select project data path
if isempty(projectPath)
    projectPath = '/Users/anouk/Documents/ExpsTablet/';
    expFolder = selectFiles([projectPath '*VMR*'],'folders');
    projectPath = [projectPath expFolder.name '/'];
end
% check if we are at the right level
while ~exist([projectPath '/1_RawData/'],'dir');
    expFolder = selectFiles(projectPath,'folders');
    projectPath = [projectPath expFolder.name '/'];
end
expName = getFolderName(projectPath); % experiment name
% define path for input and output data
dataPath = [projectPath '/1_RawData/'];
saveToPath = [projectPath '/2_ProcessedData/'];
if ~exist(saveToPath,'dir')
    mkdir(saveToPath)
end
yesGazeData = true; % default: gaze data is collected

% select subjects
cd(dataPath)
subjFolders = selectFiles('*','folders');
nSubj = length(subjFolders);

% open figure
if plotTrials
    fig1 = scaledFigure(2,1);
end

%% Loop over subjects

for s = 1 : nSubj
    
    % create structs for saving
    Exp=[];
    D=[];
    
    % get data files
    fprintf('\nLoading %s ...',subjFolders(s).name)
    subjDir = [dataPath '/' subjFolders(s).name '/'];
    controlFiles = dir([subjDir '*.txt']);
    Exp.controlFiles = {controlFiles.name}';
    fprintf(['\nControl file(s): ' repmat('%s,',1,length(controlFiles)-1) '%s\n'],Exp.controlFiles{:})
    dataFiles = dir([subjDir '*RC*.dat']);
    [~,order] = sort([dataFiles.datenum]);
    if any(diff(order)~=1)
        duplicateFiles = [dir([subjDir 'A001*RC*.dat']); dir([subjDir '*RC*copy.dat'])];
        disp('Check if trials are in correct order')
        beep; keyboard
        % to sort by date and time, execute: dataFiles = dataFiles(order,:)
    end
    nTrials = length(dataFiles);
    if nTrials==0
        disp('No data files found')
        beep; keyboard
    end
    
    % preallocate variables
    dateTime        = NaN(nTrials,6);
    trialNo         = NaN(nTrials,1);
    blockNo         = NaN(nTrials,1);
    trialType       = cell(nTrials,1);
    reachTarget     = NaN(nTrials,1);
    target1Angle    = NaN(nTrials,1);
    target1XY       = NaN(nTrials,2);
    target2Angle    = NaN(nTrials,1);
    target2XY       = NaN(nTrials,2);
    fixationAngle   = NaN(nTrials,1);
    fixationDistance = NaN(nTrials,1);
    fixationXY      = NaN(nTrials,2);
    cursorRotation  = NaN(nTrials,1);
    target1Rotation = NaN(nTrials,1);
    target2Rotation = NaN(nTrials,1);
    counter1Angle   = NaN(nTrials,1);
    counter2Angle   = NaN(nTrials,1);
    targetDistance  = NaN(nTrials,1);
    targetDistance_pen = NaN(nTrials,1);
    previewDur      = NaN(nTrials,1);
    fixationDelay   = NaN(nTrials,1);
    maxRT           = NaN(nTrials,1);
    minMT           = NaN(nTrials,1);
    maxMT           = NaN(nTrials,1);
    cursorVisible   = NaN(nTrials,1);
    feedbackMessage = cell(nTrials,1);
    tRecStartEnd    = NaN(nTrials,2);
    iTargetGoLeaveRingEnd   = NaN(nTrials,5);
    tTargetGoLeaveRingEnd   = NaN(nTrials,5);
    report1Angle_raw        = NaN(nTrials,1);
    report2Angle_raw        = NaN(nTrials,1);
    initAngle_cursor_raw    = NaN(nTrials,1);
    initAngle_pen_raw       = NaN(nTrials,1);
    hitAngle_cursor_raw     = NaN(nTrials,1);
    hitAngle_pen_raw        = NaN(nTrials,1);
    hitXY_cursor_raw        = NaN(nTrials,2);
    hitXY_pen_raw           = NaN(nTrials,2);
    fixXY_iTarget           = NaN(nTrials,2);
    D.time          = cell(nTrials,1);
    D.trialState    = cell(nTrials,1);
    D.xyPen         = cell(nTrials,1);
    D.xyCursor      = cell(nTrials,1);
    D.xyGaze_raw    = cell(nTrials,1);
    D.xyGaze        = cell(nTrials,1);
    D.vGaze         = cell(nTrials,1);
    D.dirGaze       = cell(nTrials,1);
    D.blink         = cell(nTrials,1);
    
    %% Loop over trials
    
    for t = 1 : length(dataFiles)
        
        % read in data to T structure
        if isempty(dataFiles(t).name)
            feedbackMessage{t} = 'NoDataFile';
            trialType{t} = 'NoDataFile';
            fprintf('File %d does not exist',t); continue
        end
        T = ReadFlanData([subjDir dataFiles(t).name]);
        
        %% Get landmark properties (if landmarks were present)
        % and shorten header for faster processing
        
        % for first trial of session:
        % scan header for landmark presence and get locations of landmark circles
        if t==1
            nHeaderLines = getHeaderValue(T.textdata,'HeaderLines');
            firstLandmarkLine = [];
            rectangleLine = [];
            landmarkAngle = NaN;
            landmarkRadius = NaN;
            % find first line with landmark object ('Perl' in data file)
            lineNo = 0;
            k = [];
            while isempty(k) && lineNo<=nHeaderLines
                lineNo = lineNo+1;
                k = strfind(T.textdata{lineNo},'Perl');
            end
            % if landmark is found, get properties
            if ~isempty(k)
                firstLandmarkLine = lineNo;
                i = find(isstrprop(T.textdata{firstLandmarkLine},'digit'),1); % object number
                landmarkNo = str2double(T.textdata{firstLandmarkLine}(i));
                landmarkRadius = getHeaderValue(T.textdata(firstLandmarkLine:end,1),['Object' num2str(landmarkNo) 'Radius']); % radius
                i = 1;
                for lineNo = lineNo : length(T.textdata)-1 % angles of all landmarks
                    k = strfind(T.textdata{lineNo},'Perl');
                    if ~isempty(k)
                        landmarkLine = lineNo;
                        landmarkAngle(i,1) = getHeaderValue(T.textdata(landmarkLine:end,1),['Object' num2str(landmarkNo) 'Angle']);
                        landmarkNo = landmarkNo+1;
                        i = i+1;
                    end
                end
                % find last line with landmark or text rectangle (i.e., landmark number on subject display) object
                lineNo = landmarkLine;
                for lineNo = lineNo : length(T.textdata)-1
                    k = [strfind(T.textdata{lineNo},'Perl') strfind(T.textdata{lineNo},'TextRectangle')];
                    if ~isempty(k)
                        rectangleLine = lineNo;
                    end
                end
            end
        end
        
        % shorten header to speed up code
        if ~isempty(firstLandmarkLine)
            textdata = T.textdata([1:firstLandmarkLine-1,rectangleLine(end):end-1],1);
        else
            textdata = T.textdata(1:length(T.textdata)-1,1);
        end
        textdata{1} = ['Headerlines ' num2str(length(textdata))];
        
        %% Get general experiment information
        
        % get general info from datafile of first trial
        if t==1
            
            Exp.expName = expName;
            Exp = getTabletExpInfo(Exp,textdata);
            tabletCenterXY = Exp.setup.tabletCenterXY;
            Exp.subjFolder = subjFolders(s).name;
            Exp.rawdataFiles = rmfield(dataFiles,{'bytes','isdir','datenum'});
            
            % experiment can have 1 or 2 targets, get numbers of target objects
            % search for objectname=Target1
            targetText = [];
            lineNo = 0;
            while isempty(targetText)
                lineNo = lineNo+1;
                targetText = regexp(textdata{lineNo},'Object\w*Name\sTarget1','match');
            end
            i = find(isstrprop(targetText{:},'digit'),1); % object number
            target1ObjectNo_str = targetText{:}(i);
            % search for objectname=Target2
            targetText = [];
            while isempty(targetText)
                lineNo = lineNo+1;
                targetText = regexp(textdata{lineNo},'Object\w*Name\sTarget2','match');
            end
            i = find(isstrprop(targetText{:},'digit'),1); % object number
            target2ObjectNo_str = targetText{:}(i);
            % search for objectname=fixation
            targetText = [];
            fixationExist = false;
            while isempty(targetText) && lineNo<length(textdata)
                lineNo = lineNo+1;
                targetText = regexp(textdata{lineNo},'Object\w*Name\sFixation','match');
            end
            if ~isempty(targetText)
                i = find(isstrprop(targetText{:},'digit'),1); % object number
                fixationObjectNo_str = targetText{:}(i);
                fixationExist = true;
            end
            
            % get method for determining movement endpoint
            % default is based on movement distance, some projects use movement velocity
            endpointBasedOnVelocity = nanmax([0 getHeaderValue(textdata,'ReachDetectVelocityEnable')]);
            
        end
        
        %% Get trial and stimulus information
        
        % get trial info
        dateTimeStr = getHeaderValue(textdata,'Created');
        dateTime(t,:) = datevec(dateTimeStr(2:end-1),'ddd mmm dd HH:MM:SS yyyy');
        blockNo(t) = getHeaderValue(textdata,'BlockNumber');
        trialType{t} = getHeaderValue(textdata,'TrialType');
        reportTrial = strcmp(trialType{t},'ReachDM_G') | strcmp(trialType{t},'ReachDM_I') | ...
            strcmp(trialType{t},'ReachDM_K') | strcmp(trialType{t},'ReachDM_L');
        
        % make block numbers consistent across experiments
        if strcmpi(expName,'Exp3') && blockNo(t)>1
            blockNo(t) = blockNo(t)+1;
        end
        if strcmpi(expName,'Exp3_VMR45_WR') && blockNo(t)==6
            blockNo(t) = 3;
        end
        
        % make sure trial number is correct even if block is restarted
        if t==1 || blockNo(t-1)~=blockNo(t)
            trialNo(t) = 1;
        else
            trialNo(t) = trialNo(t-1)+1;
        end
        
        % get visuomotor rotation(s)
        cursorRotation(t) = getHeaderValue(textdata,'CursorRotation');
        target1Rotation(t) = getHeaderValue(textdata,'Target1Rotation');
        target2Rotation(t) = getHeaderValue(textdata,'Target2Rotation');
        
        % get visual target angle (0 to 360 deg) and position (mm)
        reachTargetNo_str = num2str(getHeaderValue(textdata,'TargetObject'));
        if strcmp(target1ObjectNo_str,reachTargetNo_str)
            reachTarget(t) = 1;
        elseif strcmp(target2ObjectNo_str,reachTargetNo_str)
            reachTarget(t) = 2;
        end
        % target 1
        target1Visible = getHeaderValue(textdata,['Object' target1ObjectNo_str 'VisualEnable']);
        if target1Visible
            target1Angle(t) = getHeaderValue(textdata,['Object' target1ObjectNo_str 'Angle']);
            if target1Angle(t)<0
                target1Angle(t) = target1Angle(t)+360;
            end
            target1XY(t,1) = getHeaderValue(textdata,['Object' target1ObjectNo_str 'PosX'])-tabletCenterXY(1);
            target1XY(t,2) = (getHeaderValue(textdata,['Object' target1ObjectNo_str 'PosY'])-tabletCenterXY(2))*(-1); % default: down is positive
            % compute angle that would counteract the rotation
            if ~isnan(target1Rotation(t)) && target1Rotation(t)~=0 % 2-target experiment
                counter1Angle(t) = target1Angle(t) - target1Rotation(t);
            elseif cursorRotation(t)~=0 % 1 target experiments
                counter1Angle(t) = target1Angle(t) - cursorRotation(t);
            end
        end
        % target 2
        target2Visible = getHeaderValue(textdata,['Object' target2ObjectNo_str 'VisualEnable']);
        if target2Visible
            target2Angle(t) = getHeaderValue(textdata,['Object' target2ObjectNo_str 'Angle']);
            if target2Angle(t)<0
                target2Angle(t) = target2Angle(t)+360;
            end
            target2XY(t,1) = getHeaderValue(textdata,['Object' target2ObjectNo_str 'PosX'])-tabletCenterXY(1);
            target2XY(t,2) = (getHeaderValue(textdata,['Object' target2ObjectNo_str 'PosY'])-tabletCenterXY(2))*(-1); % default: down is positive
            % compute angle that would counteract the rotation
            if target2Rotation(t)~=0
                counter2Angle(t) = target2Angle(t)-target2Rotation(t);
            end
        end
        targetDistance(t) = getHeaderValue(textdata,['Object' reachTargetNo_str 'Distance']);
        targetDistance_pen(t) = targetDistance(t)/Exp.stim.cursorScalingXY(1);
        % fixation
        if fixationExist
            fixationVisible = getHeaderValue(textdata,['Object' fixationObjectNo_str 'VisualEnable']);
            if fixationVisible
                fixationAngle(t) = getHeaderValue(textdata,['Object' fixationObjectNo_str 'Angle']);
                if fixationAngle(t)<0
                    fixationAngle(t) = fixationAngle(t)+360;
                end
                fixationDistance(t) = getHeaderValue(textdata,['Object' fixationObjectNo_str 'Distance']);
                fixationXY(t,1) = getHeaderValue(textdata,['Object' fixationObjectNo_str 'PosX'])-tabletCenterXY(1);
                fixationXY(t,2) = (getHeaderValue(textdata,['Object' fixationObjectNo_str 'PosY'])-tabletCenterXY(2))*(-1); % default: down is positive
            end
        end
        
        % get trial timing
        previewDur(t) = getHeaderValue(textdata,'GoDelay');
        fixationDelay(t) = getHeaderValue(textdata,'FixationDelay');
        maxRT(t) = getHeaderValue(textdata,'MaxRT');
        minMT(t) = getHeaderValue(textdata,'MinReachDuration');
        maxMT(t) = getHeaderValue(textdata,'MaxReachDuration');
        
        % get trial feedback
        cursorVisible(t) = getHeaderValue(textdata,'Object2VisualEnable');
        feedbackMessage{t} = getHeaderValue(textdata,'FeedbackState');
        col = getColumnIndex(T.colheaders,'TrialState');
        trialState = T.data(:,col);
        if any(trialState==14) % trial was aborted
            if any(trialState==13) % reach complete before timeout
                feedbackMessage{t} = 'Good';
            elseif ~any(trialState==13) % reach not complete before timeout
                feedbackMessage{t} = 'TimeOut';
            end
        end
        while isnan(feedbackMessage{t})
            disp('No feedback message, check T.data'); beep; keyboard
        end
        % include all trials with correct timing (both hit and miss)
        feedbackGood = strcmpi(feedbackMessage{t},'Good') | strcmpi(feedbackMessage{t},'Short') | ...
            strcmpi(feedbackMessage{t},'Hit') | strcmpi(feedbackMessage{t},'Miss');
        
        %% Select raw data from target onset until end of trial
        
        % get time at start and end end of recording
        tRecStartEnd(t,1) = T.data(1,1);
        tRecStartEnd(t,2) = T.data(end,1);
        
        % remove data recorded before any target is visible (align to display delay)
        iStart = find(trialState==11,1); % start of display delay
        if isempty(iStart); iStart=1; end
        data = T.data(iStart:end,:);
        
        % remove duplicate samples if present
        col = getColumnIndex(T.colheaders,'Time');
        time_tr = data(:,col);
        notNew = [0; diff(time_tr)==0];
        data = data(~notNew,:);
        
        %% Get raw data: time, trial state, pen, cursor, gaze, and participant's report
        
        % preallocate variables
        nSamples    = size(data,1);
        xyPen       = NaN(nSamples,2);
        xyCursor    = NaN(nSamples,2);
        xyGaze_raw  = NaN(nSamples,2);
        
        % get time
        col = getColumnIndex(T.colheaders,'Time');
        time_tr = data(:,col);
        
        % get trial state names and durations
        col = getColumnIndex(T.colheaders,'TrialState');
        trialState = data(:,col);
        trialStateNames = Exp.timing.trialStateNames;
        trialStateDur = NaN(length(trialStateNames),1);
        for i = 1 : length(trialStateNames)
            iStateOn = find(trialState==i,1);
            iStateOff = find(trialState==i,1,'last');
            if ~isempty(iStateOn)
                trialStateDur(i) = time_tr(iStateOff) - time_tr(iStateOn) + 1/Exp.setup.fs;
            end
        end
        
        % get timing of events (from trial states)
        iTarget = min([find(trialState==6,1); find(trialState==9,1); find(trialState==7,1)]); % display scene or report or fill target
        iGo = find(trialState==7,1);        % fill target is go cue
        iLeave = find(trialState==12,1);    % detection of movement start (online)
        iComplete = find(trialState==13,1); % reach completed (online)
        iEndOfTrial = length(trialState);   % end of recording
        if isempty(iTarget) && feedbackGood
            fprintf('No ShowScene state in block %d trial %d ',blockNo(t),trialNo(t))
            iTarget = 1;
        end
        if isempty(iGo) && feedbackGood
            fprintf('No FillTarget trial state in block %d trial %d',blockNo(t),trialNo(t))
            iGo = iLeave;
        end
        
        % get pen on tablet position
        col = getColumnIndex(T.colheaders,'TabletX');
        xyPen(:,1) = data(:,col);
        col = getColumnIndex(T.colheaders,'TabletY');
        xyPen(:,2) = data(:,col)*(-1);
        
        % get cursor on screen position
        col = getColumnIndex(T.colheaders,'CursorPosX');
        xyCursor(:,1) = data(:,col)-tabletCenterXY(1);
        col = getColumnIndex(T.colheaders,'CursorPosY');
        xyCursor(:,2) = (data(:,col)-tabletCenterXY(2))*(-1);
        
        % on fMRI-compatible tablet, sampling issues can occur when the
        % finger doesn't sufficiently press down on the tablet,
        % classify as sampling issue if no samples were collected during reach
        xyCursor_reach = xyCursor(trialState==12,:);
        if feedbackGood && isempty(xyCursor_reach)
            feedbackMessage{t} = 'SamplingIssue';
            feedbackGood = false;
        end
        
        % get gaze position if gaze data was collected
        if yesGazeData
            col = getColumnIndex(T.colheaders,'GazeX');
            if col>0
                xyGaze_raw(:,1) = data(:,col)-tabletCenterXY(1);
                col = getColumnIndex(T.colheaders,'GazeY');
                xyGaze_raw(:,2) = (data(:,col)-tabletCenterXY(2))*(-1); % default: down is positive
            else
                yesGazeData = false; % no gaze data was collected
            end
        end
        
        % in report trials, get participant's verbal or dial report
        if reportTrial
            input1Value = getHeaderValue(textdata,'Input1Value');
            input1Button = getHeaderValue(textdata,'Input1Button');
            if input1Button==0 % check if button was pressed (fMRI setup)
                input1Value = NaN;
            end
            input2Value = getHeaderValue(textdata,'Input2Value'); % 2-target experiment
            % report 1
            if ~isnan(input1Value)
                if strcmp(trialType{t},'ReachDM_K') || strcmp(trialType{t},'ReachDM_L') % dial trials
                    report1Angle_raw(t) = input1Value;
                elseif isnumeric(input1Value) && input1Value<=length(landmarkAngle)
                    % convert number to angle (0 to 360 deg)
                    if ~isnan(target1Angle(t)) % report corresponds to target1
                        report1Angle_raw(t) = landmarkAngle(input1Value);
                    elseif ~isnan(target2Angle(t)) % report corresponds to target2
                        report2Angle_raw(t) = landmarkAngle(input1Value);
                    end
                else
                    fprintf('No report1 in block %d trial %d, value is %s',blockNo(t),trialNo(t),input1Value)
                end
            end
            % report 2
            if ~isnan(input2Value)
                if isnumeric(input2Value) && input2Value<=length(landmarkAngle)
                    % convert number to angle (0 to 360 deg)
                    report2Angle_raw(t) = landmarkAngle(input2Value);
                else
                    fprintf('No report1 in block %d trial %d, value is %s',blockNo(t),trialNo(t),input1Value)
                end
            end
            
        end % if reportTrial
        
        % compute cursor distance and direction during reach movement
        dCursor = sqrt(xyCursor(:,1).^2+xyCursor(:,2).^2); % distance from start
        iHalfway = find(dCursor>=(targetDistance(t)/2),1);
        iRing = find(dCursor>=targetDistance(t),1); % first sample after crossing the ring
        if endpointBasedOnVelocity
            iRing = iComplete;
        end
        
        % keep trials that were categorized as 'Short' online, but cursor
        % reached the ring, if movement time was within the limit
        if strcmpi(feedbackMessage{t},'Short')
            if isempty(iRing) % cursor didn't reach the target
                feedbackGood = false;
            else % cursor reached the target
                MT = time_tr(iRing) - time_tr(iLeave);
                if MT<maxMT(t)
                    feedbackMessage{t} = 'Good';
                else
                    feedbackGood = false;
                end
            end
        end
        
        % if timing is incorrect, move on to next trial
        if ~feedbackGood
            continue
        end
        
        %% Compute cursor direction and hit angle
        
        iRing0 = find(dCursor==dCursor(iRing-1),1); % last new sample before crossing the ring
        
        % compute initial movement angle (when cursor is halfway the target distance)
        if ~isempty(iHalfway)
            iHalfway0 = find(dCursor==dCursor(iHalfway-1),1);
            
            % x and y position of cursor
            time_interp = fliplr(time_tr(iHalfway):-0.001:time_tr(iHalfway0));
            xCross = interp1(time_tr([iHalfway0,iHalfway]),xyCursor([iHalfway0,iHalfway],1),time_interp);
            yCross = interp1(time_tr([iHalfway0,iHalfway]),xyCursor([iHalfway0,iHalfway],2),time_interp);
            dCross = sqrt(xCross.^2+yCross.^2);
            xy = [xCross(find(dCross>=(targetDistance(t)/2),1)) yCross(find(dCross>=(targetDistance(t)/2),1))];
            initAngle_cursor_raw(t) = xy2compassAngle(xy(1),xy(2),360);
            
            % x and y position of pen
            time_interp = fliplr(time_tr(iHalfway):-0.001:time_tr(iHalfway0));
            xCross = interp1(time_tr([iHalfway0,iHalfway]),xyPen([iHalfway0,iHalfway],1),time_interp);
            yCross = interp1(time_tr([iHalfway0,iHalfway]),xyPen([iHalfway0,iHalfway],2),time_interp);
            dCross = sqrt(xCross.^2+yCross.^2);
            xy = [xCross(find(dCross>=(targetDistance_pen(t)/2),1)) yCross(find(dCross>=(targetDistance_pen(t)/2),1))];
            initAngle_pen_raw(t) = xy2compassAngle(xy(1),xy(2),360);
        else
            fprintf('Block %d trial %d - Cursor did not reach 50% of target distance',blockNo(t),trialNo(t));
            keyboard
            feedbackMessage{t} = 'Short';
        end
        
        % fMRI: check distance of last sample before crossing the ring
        % classify as sampling issue if distance <50%
        if strcmpi(Exp.setup.location,'fMRI') && dCursor(iRing0)<Exp.timing.leaveStartDistance;
            iRing=[]; iRing0 = [];
            feedbackMessage{t} = 'SamplingIssue';
            hitXY_cursor_raw(t,:) = [NaN NaN];
            hitXY_pen_raw(t,:) = [NaN NaN];
        end
        
        % compute hit position and angle (when cursor crosses the ring)
        if ~endpointBasedOnVelocity
            % get x and y position of cursor when crossing the ring
            if strcmpi(Exp.setup.location,'fMRI') % fMRI: use last sample before crossing the ring
                hitXY_cursor_raw(t,:) = [xyCursor(iRing0,1) xyCursor(iRing0,2)];
            else % interpolate cursor position when crossing the ring
                time_interp = fliplr(time_tr(iRing):-0.001:time_tr(iRing0));
                xCross = interp1(time_tr([iRing0,iRing]),xyCursor([iRing0,iRing],1),time_interp);
                yCross = interp1(time_tr([iRing0,iRing]),xyCursor([iRing0,iRing],2),time_interp);
                dCross = sqrt(xCross.^2+yCross.^2);
                hitXY_cursor_raw(t,:) = [xCross(find(dCross>=targetDistance(t),1)) yCross(find(dCross>=targetDistance(t),1))];
            end
            
            % get x and y position of pen when crossing the ring
            if strcmpi(Exp.setup.location,'fMRI') % fMRI: use last sample before crossing the ring
                hitXY_pen_raw(t,:) = [xyPen(iRing0,1) xyPen(iRing0,2)];
            else % interpolate pen position when crossing the ring
                time_interp = fliplr(time_tr(iRing):-0.001:time_tr(iRing0));
                xCross = interp1(time_tr([iRing0,iRing]),xyPen([iRing0,iRing],1),time_interp);
                yCross = interp1(time_tr([iRing0,iRing]),xyPen([iRing0,iRing],2),time_interp);
                dCross = sqrt(xCross.^2+yCross.^2);
                hitXY_pen_raw(t,:) = [xCross(find(dCross>=targetDistance_pen(t),1)) yCross(find(dCross>=targetDistance_pen(t),1))];
            end
            
        elseif endpointBasedOnVelocity==1
            % get x and y position of cursor and pen
            hitXY_cursor_raw(t,:) = xyCursor(iRing,:);
            hitXY_pen_raw(t,:) = xyPen(iRing,:);
        end
        
        % get angle at hit
        hitAngle_cursor_raw(t) = xy2compassAngle(hitXY_cursor_raw(t,1),hitXY_cursor_raw(t,2),360);
        hitAngle_pen_raw(t) = xy2compassAngle(hitXY_pen_raw(t,1),hitXY_pen_raw(t,2),360);
        
        % save timing
        iTargetGoLeaveRingEnd(t,:) = [iTarget iGo iLeave iRing iEndOfTrial];
        tTargetGoLeaveRingEnd(t,:) = time_tr(iTargetGoLeaveRingEnd(t,:))-time_tr(1);
        
        %% Gaze data: remove blinks, low-pass filter and compute velocity
        
        if yesGazeData
            [xyGaze,vxyGaze,blink] = filterGazeData(xyGaze_raw,Exp);
            vGaze = sqrt(vxyGaze(:,1).^2 + vxyGaze(:,2).^2); % resultant velocity
            
            % compute gaze position at the time of target presentation
            fixXY = nanmean(xyGaze(iTarget:iTarget+19,:));
            if all(abs(fixXY)<(0.5*targetDistance(t)))
                fixXY_iTarget(t,:) = fixXY;
            else
                fixXY_iTarget(t,:) = [NaN NaN];
            end
            
            % compute gaze direction (deg) and distance
            dirGaze = xy2compassAngle(xyGaze(:,1),xyGaze(:,2),360);
            distGaze = sqrt(xyGaze(:,1).^2 + xyGaze(:,2).^2);
            
        else
            xyGaze = xyGaze_raw;
            vGaze = [];
            dirGaze = [];
            distGaze = [];
            blink = [];
        end
        
        %% Store time, cursor data, and gaze data in struct
        
        D.time{t,1}         = time_tr;
        D.trialState{t,1}   = trialState;
        D.xyCursor{t,1}     = xyCursor;
        D.xyPen{t,1}        = xyPen;
        D.xyGaze_raw{t,1}   = xyGaze_raw;
        D.xyGaze{t,1}       = xyGaze;
        D.vGaze{t,1}        = vGaze;
        D.dirGaze{t,1}      = dirGaze;
        D.blink{t,1}        = blink;
        
        %% Plot raw data of trials with good timing
        
        if plotTrials && cursorRotation(t)~=0
            d = targetDistance(t);
            MT = tTargetGoLeaveRingEnd(t,4)-tTargetGoLeaveRingEnd(t,3); % movement time
            % plot data in xy coordinates
            figure(fig1); %clf;
            subplot(1,2,1); cla
            int = iTarget:nanmin([iRing nSamples]); % interval from target onset to hit (or end of trial)
            plot(d*sind(landmarkAngle),d*cosd(landmarkAngle),'o','color',[0.5 0.5 0.5]); % plot ring of pearls
            hold on
            plot(landmarkRadius*sind(0:360),landmarkRadius*cosd(0:360),'k'); % plot single ring
            pg = plot(xyGaze(int,1),xyGaze(int,2),'.-','color',[0 0.5 0]); % plot gaze
            if cursorVisible(t)
                pc = plot(xyCursor(int,1),xyCursor(int,2),'.-','color','c'); % plot cursor
            else
                pc = plot(xyCursor(int,1),xyCursor(int,2),'.-','color','y');
            end
            pt1 = plot(target1XY(t,1),target1XY(t,2),'ro','markerfacecolor','r','markersize',8); % plot target 1
            plot([0 target1XY(t,1)],[0 target1XY(t,2)],'r');
            pt2 = plot(target2XY(t,1),target2XY(t,2),'bo','markerfacecolor','b','markersize',8); % plot target 2
            plot([0 target2XY(t,1)],[0 target2XY(t,2)],'b')
            pf = plot(fixationXY(t,1),fixationXY(t,2),'k+');                                    % plot fixation
            plot(hitXY_cursor_raw(t,1),hitXY_cursor_raw(t,2),'ko')                                      % plot hit
            pa1 = plot([0 d*sind(counter1Angle(t))],[0 d*cosd(counter1Angle(t))],'r--');        % plot aimpoint 1
            pa2 = plot([0 d*sind(counter2Angle(t))],[0 d*cosd(counter2Angle(t))],'b--');        % plot aimpoint 2
            pr1 = plot([0 d*sind(report1Angle_raw(t))],[0 d*cosd(report1Angle_raw(t))],'r:');   % plot report 1
            pr2 = plot([0 d*sind(report2Angle_raw(t))],[0 d*cosd(report2Angle_raw(t))],'b:');   % plot report 2
            % axes
            if any(~isnan(target2XY))
                legend([pc pg pt1 pt2 pa1 pa2 pr1 pr2],...
                    {'cursor','gaze','target1','target2','aimpoint1','aimpoint2','report1','report2'});
            else
                legend([pc pg pt1 pa1 pr1],...
                    {'cursor','gaze','target','aimpoint','report'});
            end
            legend('boxoff')
            axis(d*[-1.15 1.35 -1.15 1.35]); axis square
            text(d*-1.05,d*1.2,['MT = ' num2str(round(MT*1000)) ' ms'])
            xlabel('X (mm)'); ylabel('Y (mm)')
            title(['Block ' num2str(blockNo(t)) ' trial ' num2str(trialNo(t))])
            hold off
            
            % plot gaze direction over time
            if yesGazeData
                subplot(1,2,2); cla reset
                angles = [target1Angle(t) target2Angle(t) counter1Angle(t) counter2Angle(t) report1Angle_raw(t) report2Angle_raw(t) hitAngle_cursor_raw(t)-cursorRotation(t)];
                angles(2,:) = angles;
                angles(2,angles(1,:)<=0) = angles(1,angles(1,:)<=0)+360;
                xlim([time_tr(iTarget) time_tr(iRing)])
                horline(angles(:,1),'r')    % target 1
                horline(angles(:,2),'b')    % target 2
                horline(angles(:,3),'r--')  % counter 1
                horline(angles(:,4),'b--')  % counter 2
                horline(angles(:,5),'r:')   % report 1
                horline(angles(:,6),'b:')   % report 2
                horline(angles(:,7),'c-.')  % hand direction at hit
                goodDistance = distGaze>(d*0.5) & distGaze<(d*1.5);
                hold on; plot(time_tr(goodDistance),dirGaze(goodDistance),'o','color',[0 0.5 0]) % gaze
                hold off
                ylim([-5 365])
                if ~isempty(iGo)
                    vertline(time_tr(iGo))
                    text(time_tr(iGo),350,'GO')
                    vertline(time_tr(iLeave))
                    text(time_tr(iLeave),350,'reach')
                end
                set(gca,'Ytick',0:45:360)
                title('Gaze direction (at 50-150% of target distance)')
                xlabel('Time since target onset (s)')
                ylabel('Direction (deg)')
            end
            waitforbuttonpress % click/press key to continue
        end
        
    end % end of loop over trials
    
    %% All trials: compute cursor angle (=error) and report angle relative to target angle
    
    % compute directional error: cursor angle relative to target angle (deg)
    initAngle_cursor = NaN(nTrials,2);
    initAngle_cursor(reachTarget==1,1) = initAngle_cursor_raw(reachTarget==1) - target1Angle(reachTarget==1);
    initAngle_cursor(reachTarget==2,2) = initAngle_cursor_raw(reachTarget==2) - target2Angle(reachTarget==2);
    initAngle_cursor(initAngle_cursor>180) = initAngle_cursor(initAngle_cursor>180)-360;
    initAngle_cursor(initAngle_cursor<-180) = initAngle_cursor(initAngle_cursor<-180)+360;
    hitAngle_cursor = NaN(nTrials,2);
    hitAngle_cursor(reachTarget==1,1) = hitAngle_cursor_raw(reachTarget==1) - target1Angle(reachTarget==1);
    hitAngle_cursor(reachTarget==2,2) = hitAngle_cursor_raw(reachTarget==2) - target2Angle(reachTarget==2);
    hitAngle_cursor(hitAngle_cursor>180) = hitAngle_cursor(hitAngle_cursor>180)-360;
    hitAngle_cursor(hitAngle_cursor<-180) = hitAngle_cursor(hitAngle_cursor<-180)+360;
    
    % compute directional error: pen angle relative to target angle (deg)
    initAngle_pen = NaN(nTrials,2);
    initAngle_pen(reachTarget==1,1) = initAngle_pen_raw(reachTarget==1) - target1Angle(reachTarget==1);
    initAngle_pen(reachTarget==2,2) = initAngle_pen_raw(reachTarget==2) - target2Angle(reachTarget==2);
    initAngle_pen(initAngle_pen>180) = initAngle_pen(initAngle_pen>180)-360;
    initAngle_pen(initAngle_pen<-180) = initAngle_pen(initAngle_pen<-180)+360;
    hitAngle_pen = NaN(nTrials,2);
    hitAngle_pen(reachTarget==1,1) = hitAngle_pen_raw(reachTarget==1) - target1Angle(reachTarget==1);
    hitAngle_pen(reachTarget==2,2) = hitAngle_pen_raw(reachTarget==2) - target2Angle(reachTarget==2);
    hitAngle_pen(hitAngle_pen>180) = hitAngle_pen(hitAngle_pen>180)-360;
    hitAngle_pen(hitAngle_pen<-180) = hitAngle_pen(hitAngle_pen<-180)+360;
    
    % compute hand angle relative to target angle
    initAngle_hand = initAngle_cursor - [cursorRotation cursorRotation];
    initAngle_hand(initAngle_hand<-180) = initAngle_hand(initAngle_hand<-180)+360;
    hitAngle_hand = hitAngle_cursor - [cursorRotation cursorRotation];
    if str2double(expName(isstrprop(expName,'digit')))==4
        hitAngle_hand = hitAngle_pen;
    end
    hitAngle_hand(hitAngle_hand<-180) = hitAngle_hand(hitAngle_hand<-180)+360;
    
    % compute explicit component: report angle relative to target angle
    explicit1Angle = report1Angle_raw - target1Angle;
    explicit1Angle(explicit1Angle>180) = explicit1Angle(explicit1Angle>180)-360;
    explicit1Angle(explicit1Angle<-180) = explicit1Angle(explicit1Angle<-180)+360;
    explicit2Angle = report2Angle_raw - target2Angle;
    explicit2Angle(explicit2Angle>180) = explicit2Angle(explicit2Angle>180)-360;
    explicit2Angle(explicit2Angle<-180) = explicit2Angle(explicit2Angle<-180)+360;
    
    % implicit component can be calculated after binning the trials
    
    %% Print number of good trials and error trials in command window
    
    nGood = sum(ismember(feedbackMessage,'Good'));
    nEarly = sum(ismember(feedbackMessage,'Early'));
    nLate = sum(ismember(feedbackMessage,'Late'));
    nShort = sum(ismember(feedbackMessage,'Short') | ismember(feedbackMessage,'ReachTimeOut'));
    nTimeOut = sum(ismember(feedbackMessage,'TimeOut'));
    nSamplingIssues = sum(ismember(feedbackMessage,'SamplingIssue'));
    fprintf('\n%d / %d trials with good timing',nGood,nTrials)
    fprintf('\n%d trials early',nEarly)
    fprintf('\n%d trials late',nLate)
    fprintf('\n%d trials short or slow',nShort)
    fprintf('\n%d trials with timeout',nTimeOut)
    fprintf('\n%d trials with sampling issues',nSamplingIssues)
    
    nReport = sum(~isnan(explicit1Angle));
    nReportTrials = sum(ismember(trialType,'ReachDM_G') | ismember(trialType,'ReachDM_K') | ismember(trialType,'ReachDM_L'));
    fprintf('\n%d / %d trials with report\n',nReport,nReportTrials);
    
    %% Create structs for saving
    
    Exp.timing.previewDur = previewDur;
    Exp.timing.maxRT    = maxRT;
    Exp.timing.minMT    = minMT;
    Exp.timing.maxMT    = maxMT;
    Exp.processingDate  = date;
    Exp.processingCode  = 'processVMRtabletData.m by Anouk de Brouwer';
    Exp.blockNo         = blockNo;
    Exp.trialNo         = trialNo;
    Exp.trialType       = trialType;
    
    % stimulus information
    Exp.reachTarget     = reachTarget;
    Exp.targetDistance  = targetDistance;
    Exp.target1XY       = target1XY;
    Exp.target1Angle    = target1Angle;
    Exp.target2XY       = target2XY;
    Exp.target2Angle    = target2Angle;
    Exp.fixationAngle   = fixationAngle;
    Exp.fixationDistance = fixationDistance;
    Exp.fixationXY      = fixationXY;
    Exp.cursorVisible   = cursorVisible;
    Exp.cursorRotation  = cursorRotation;
    Exp.target1Rotation = target1Rotation;
    Exp.target2Rotation = target2Rotation;
    Exp.counter1Angle   = counter1Angle;
    Exp.counter2Angle   = counter2Angle;
    Exp.pearlNumAngle   = [(1:length(landmarkAngle))' landmarkAngle];
    
    % data: feedback and timing
    % already defined in trial loop:
    % D.time{t,1} = time_tr;
    % D.trialState{t,1} = trialState;
    D.feedbackMessage = feedbackMessage;
    D.iTargetGoLeaveRingEnd = iTargetGoLeaveRingEnd;
    D.tTargetGoLeaveRingEnd = tTargetGoLeaveRingEnd;
    D.fixXY_iTarget = fixXY_iTarget;
    
    % data: hand and gaze
    % already defined in trial loop:
    % D.xyPen{t,1}      = xyPen;
    % D.xyCursor{t,1}   = xyCursor;
    % D.xyGaze_raw{t,1} = xyGaze_raw;
    % D.xyGaze{t,1}     = xyGaze;
    % D.vGaze{t,1}      = vGaze;
    % D.dirGaze{t,1}    = dirGaze;
    % D.blink{t,1}      = blink;
    D.explicit1Angle    = explicit1Angle;
    D.explicit2Angle    = explicit2Angle;
    D.initAngle_hand    = initAngle_hand;
    D.initAngle_cursor  = initAngle_cursor;
    D.hitAngle_hand     = hitAngle_hand;
    D.hitAngle_cursor   = hitAngle_cursor;
    D.hitXY_hand_raw    = hitXY_pen_raw;
    D.hitXY_cursor_raw  = hitXY_cursor_raw;
    D.description = {'All data is aligned to start of displayDelay,';...
        'occurring after the hand is at start and before the stimuli are presented.';...
        'XY positions are in mm, directions are in degrees.'};
    
    %% Save
    
    % check if file exists
    fileName = [subjFolders(s).name '.mat'];
    if exist([saveToPath fileName],'file') == 2 % check if file does not exist yet
        disp(['A file named ' fileName ' already exists.'])
        overwrite = input('Do you want to overwrite it? Yes(1) or no(0): ');
    else
        overwrite = 1;
    end
    
    % save file
    if overwrite == 1
        save([saveToPath fileName],'Exp','D');
        disp(['Saved ' saveToPath fileName])
    else
        disp('Data has not been saved')
    end
    
    if plotTrials==1
        disp('Click continue to go to next subject')
        beep; keyboard
    end
    
end % end of loop over subjects

% save copy of Matlab code used to process data
mFilePath = mfilename('fullpath');
saveCopyOfCode(mFilePath,saveToPath)

end