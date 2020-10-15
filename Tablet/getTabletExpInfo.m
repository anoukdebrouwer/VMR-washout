function [Exp] = getTabletExpInfo(Exp,textdata)
% getTabletExpInfo  Get tablet experiment information from header of datafile
% 
% [Exp] = getTabletExpInfo(Exp,textdata) gets general information about 
% subject, task program, stimuli, timing and setup from the header of the 
% raw datafile and adds it to the struct Exp.

% MIT License
% Copyright (c) 2020 Anouk de Brouwer

% subject
subjStr = getHeaderValue(textdata,'Subject');
Exp.subjName = subjStr(isstrprop(subjStr,'alphanum'));

% task program
Exp.programName = getHeaderValue(textdata,'Program');
Exp.expDateAndTime = getHeaderValue(textdata,'Created');
disp(['First trial: ' Exp.expDateAndTime(2:end-1)])

%% Stimuli

% start
stim.startXY = [getHeaderValue(textdata,'StartPosX') getHeaderValue(textdata,'StartPosY')];
objectText = [];
lineNo = 0;
while isempty(objectText)
    lineNo = lineNo+1;
    objectText = regexp(textdata{lineNo},'Object\w*Name\sStart','match');
end
i = find(isstrprop(objectText{:},'digit'),1); % object number
objectNo_str = objectText{:}(i);
stim.startRadius = getHeaderValue(textdata(lineNo:end,1),['Object' objectNo_str 'Radius']);

% cursor
objectText = [];
lineNo = 0;
while isempty(objectText)
    lineNo = lineNo+1;
    objectText = regexp(textdata{lineNo},'Object\w*Name\sCursor','match');
end
i = find(isstrprop(objectText{:},'digit'),1); % object number
objectNo_str = objectText{:}(i);
stim.cursorRadius = getHeaderValue(textdata(lineNo:end,1),['Object' objectNo_str 'Radius']);           

% scaling between pen position on tablet and cursor position on screen
stim.cursorScalingXY = [getHeaderValue(textdata,'CursorScaleX'),...
    getHeaderValue(textdata,'CursorScaleY')]; 
if isnan(stim.cursorScalingXY)
    stim.cursorScalingXY = [1 1];
end

% target
objectText = [];
lineNo = 0;
while isempty(objectText)
    lineNo = lineNo+1;
    objectText = regexp(textdata{lineNo},'Object\w*Name\sTarget','match');
end
i = find(isstrprop(objectText{:},'digit'),1); % object number
objectNo_str = objectText{:}(i);
stim.targetRadius = getHeaderValue(textdata(lineNo:end,1),['Object' objectNo_str 'Radius']);
stim.targetWidthHeight = [getHeaderValue(textdata(lineNo:end,1),['Object' objectNo_str 'SizeX']),...
    getHeaderValue(textdata(lineNo:end,1),['Object' objectNo_str 'SizeY'])];
stim.targetDistance = getHeaderValue(textdata(lineNo:end,1),['Object' objectNo_str 'Distance']);

% fixation (if present)
objectText = [];
lineNo = 0;
while isempty(objectText) && lineNo<length(textdata)
    lineNo = lineNo+1;
    objectText = regexp(textdata{lineNo},'Object\w*Name\sFixation','match');
end
if ~isempty(objectText)
    i = find(isstrprop(objectText{:},'digit'),1); % object number
    objectNo_str = objectText{:}(i);
    stim.fixationWidthHeight = [getHeaderValue(textdata(lineNo:end,1),['Object' objectNo_str 'SizeX']),...
        getHeaderValue(textdata(lineNo:end,1),['Object' objectNo_str 'SizeY'])];
end

Exp.stim = stim;

%% Trial timing

% get trial state names
i=1;
if ~isnan(getHeaderValue(textdata,'TrialState[1]'))
    while ~isnan(getHeaderValue(textdata,['TrialState[' num2str(i) ']']))
        trialStateNames{i,1} = getHeaderValue(textdata,['TrialState[' num2str(i) ']']);
        i = i+1;
    end
else
    trialStateNames = {};
end
timing.trialStateNames = trialStateNames;

% trial timing
timing.holdStartDistance = getHeaderValue(textdata,'CursorAtStartDistance');
timing.holdStartDur = getHeaderValue(textdata,'HoldStartDuration');
timing.preTargetDelay = getHeaderValue(textdata,'DisplayDelay');
timing.leaveStartDistance = getHeaderValue(textdata,'LeaveStartDistance');
timing.endPointFeedbackDur = getHeaderValue(textdata,'PostReachDelay');
timing.postReachFeedbackDur = getHeaderValue(textdata,'FeedbackDuration');
timing.maxTrialDuration = getHeaderValue(textdata,'MaxTrialDuration');
timing.note = {'Feedback duration and visibility might differ per block and trial'};

Exp.timing = timing;

%% Setup

% location and sampling
setup.name = 'Tablet';
setup.fs = getHeaderValue(textdata,'DataRecordFrequency');
if ~isempty(strfind(Exp.controlFiles{1},'fmri'))
    setup.location = 'fMRI';
else
    setup.fs_tablet = 100;
    setup.fs_Eyelink = 500;
    setup.location = 'Abramsky';
end

% tablet specs
setup.tabletSize_mm = [getHeaderValue(textdata,'TabletSizeX') getHeaderValue(textdata,'TabletSizeY')];
tabletCenterXY_mm = [getHeaderValue(textdata,'TabletCenterX') getHeaderValue(textdata,'TabletCenterY')];
if isnan(tabletCenterXY_mm(1)); tabletCenterXY_mm = [0 0]; end
setup.tabletCenterXY = tabletCenterXY_mm;

% monitor specs
setup.monitorSize_mm = [getHeaderValue(textdata,'DisplaySizeX') getHeaderValue(textdata,'DisplaySizeY')];
setup.monitorCentreXY_mm = [getHeaderValue(textdata,'DisplayCenterX') getHeaderValue(textdata,'DisplayCenterY')];

Exp.setup = setup;