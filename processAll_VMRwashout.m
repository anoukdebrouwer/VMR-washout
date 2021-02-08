% processAll_VMRwashout

% MIT License
% Copyright (c) 2020 Anouk de Brouwer

projectPath = '/Users/anouk/Documents/ExpsTablet/VMRlearning_gaze/';
saveToPath = '/Users/anouk/dropbox/Gaze_VMR_adaptation/';

%% Process VMR task data of each experiment

% define data path - choose one experiment
%expPath = [projectPath 'Exp1_reporting/'];
%expPath = [projectPath 'Exp2_noReporting/'];
expPath = [projectPath 'Exp7_WR/'];

% process raw data 
plotTrials = false; % set to true to plot each trial separately for visual inspection
processVMRtabletData(expPath,plotTrials);

% calculate and plot individual subject results
processData = true;
createPlots = true;
savePlots = true;
calcIndResultsVMRgaze(expPath,processData,createPlots,savePlots)

% combine individual subject results into a single datafile
meanOrMedian = 1;
createPlots = true;
savePlots = true;
combineIndResultsVMRgaze(expPath,meanOrMedian,createPlots,savePlots)

%% Combine the data of all experiments into a single table 

createResultsTable_VMRwashout(projectPath,saveToPath)

%% Final results

% See VMRwashout.ipynb
