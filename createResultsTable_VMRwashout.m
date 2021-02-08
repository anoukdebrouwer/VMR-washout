function createResultsTable_VMRwashout(projectPath,saveToPath)
% creatResultsTable_VMRwashout Create table with results of VMR washout project
%
% createResultsTable_VMRwashout(projectPath,saveToPath) loads the processed
% data file of each experiment in projectPath selected by the user, and
% combines them in a single table. The table is saved as a .csv file in
% saveToPath.

% MIT License
% Copyright (c) 2020 Anouk de Brouwer

if nargin==0
    projectPath = '/Users/anouk/Documents/ExpsTablet/VMRlearning_gaze/';
    saveToPath = '/Users/anouk/dropbox/Gaze_VMR_adaptation/';
end

% select experiment folders
expFolders = selectFiles([projectPath 'Exp*'],'folders');

% loop over experiment folders
resultsFiles = [];
for exp = 1 : length(expFolders)
    dataPath = [projectPath expFolders(exp).name '/3_Results/'];
    % select subject data files
    selectedFiles = selectFiles([dataPath 'Results*.mat'],'files');
    resultsFiles = [resultsFiles; strcat(dataPath,{selectedFiles.name}')];
end
nFiles = length(resultsFiles);

% loop over experiment data files
T = table;
for f = 1 : nFiles
    
    % load file
    load(resultsFiles{f})
    disp(['Loaded ' resultsFiles{f}])
    expName = getFolderName(resultsFiles{f},2);
    subj_all = subj;
    nSubj = length(subj_all);
    
    %% Create tables
    
    % concatenate experiment specs across 2 days
    nBins = size(percReportStrategy,1);
    nRows = nBins*2;
    exp = repmat({expName},nRows,1);
    day = [ones(nBins,1); ones(nBins,1)*2];
    uniBlockName = unique(blockName,'stable');
    bin_no = [];
    for i = 1 : length(uniBlockName)
        nBins_block = sum(ismember(blockName,uniBlockName(i)));
        bin_no = [bin_no 1:nBins_block];
    end
    block_name = [blockName(:); blockName(:)];
    bin_no = [bin_no(:); bin_no(:)];
    report = [report(:); report(:)];
    
    % add data of each subject to table
    for s = 1 : nSubj
        % check if subject code is unique
        if f>1
            if ismember(subj_all{s},T.subj);
                alpha = isstrprop(subj_all{s},'alpha');
                num = isstrprop(subj_all{s},'digit');
                subj_new = [subj_all{s}(alpha) num2str(str2double(subj_all{s}(num))+1)];
                disp([subj_all{s} ' was renamed to ' subj_new])
                subj_all{s} = subj_new;
            end
        end
        % reshape variables for long table format
        subj = repmat(subj_all(s),nRows,1);
        hand_angle = reshape(handAngle(:,s,1:2),[],1);
        report_angle = reshape(reportAngle(:,s,1:2),[],1);
        aimfix_angle = reshape(aimFixAngle(:,s,1:2),[],1);
        implicit_angle = reshape(implicitAngle(:,s,1:2),[],1);
        implicit_angle_from_fix = reshape(implicitAngle_fromFix(:,s,1:2),[],1);
        perc_report_strategy = reshape(percReportStrategy(:,s,1:2),[],1);
        perc_fix_strategy = reshape(percFixStrategy(:,s,1:2),[],1);
        % create table
        T_subj = table(exp,day,block_name,bin_no,subj,...
            hand_angle,report_angle,aimfix_angle,...
            implicit_angle,implicit_angle_from_fix,...
            perc_report_strategy,perc_fix_strategy);
        T = [T; T_subj];
    end
    
end

% check if file exists and save
fileName = ['table_VMRwashout_' datestr(now,'yyyymm') '.csv'];
overwrite = checkFileExistence(saveToPath,fileName);
if overwrite == 1
    writetable(T,[saveToPath fileName],'WriteRowNames',true);
    disp(['Saved ' saveToPath fileName])
else
    disp('Table has not been saved')
    keyboard
end