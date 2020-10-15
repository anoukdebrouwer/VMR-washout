function combineIndResultsVMRgaze(projectPath,meanOrMedian,createPlots,savePlots)
% combineIndResultsVMRgaze Bin individual data of VMR gaze experiments and
% combine into a group data file.
%
% combineResultsVMRgaze(projectPath,meanOrMedian) loads the individual data
% files in projectPath, calculates the mean (1; default) or median (2)
% hand angle, reported aiming angle, aim fixation angle, and implicit angle
% per bin of trials for each participant. It then saves a matrix with the
% for each of these variables with the values of all participants.
%
% combineResultsVMRgaze(projectPath,meanOrMedian,createPlots) plots binned
% angles when createPlots=TRUE and saves these plots when savePlots=true;.
%
% For 2-day experiments, it is assumed that the experimental blocks on both
% days are identical (although day 2 doesn't need to include all day-1 blocks).

% MIT License
% Copyright (c) 2020 Anouk de Brouwer

if nargin==0
    projectPath = [];
    meanOrMedian = 1;
    createPlots = true;
    savePlots = false;
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
expName = getFolderName(projectPath);
% define input and output data path
dataPath = [projectPath '/3_Results/'];
saveToPath = [projectPath '/3_Results/'];
cd(dataPath)

% select subject data files
subjFiles = selectFiles('S*.mat');
if isempty(subjFiles)
    subjFiles = selectFiles('*.mat');
end
nSubj = length(subjFiles);
subj = cell(1,nSubj);

% load experiment details
detailsFile = dir([projectPath '/ExpDetails*.mat']);
load([projectPath '/' detailsFile.name])
nBins = nTrials/nTargets;
vmr = visuomotorRotation;
dPearlAngle = abs(pearlNumAngle(2,2)-pearlNumAngle(1,2));
targetZone = [-1.5*dPearlAngle 1.5*dPearlAngle];

% open figures
fig1 = figure;
matlabcolors = get(gca,'colororder');

% preallocate
RT                      = NaN(nBins,nSubj,nDays);
sdRT                    = NaN(nBins,nSubj,nDays);
cursorAngle             = NaN(nBins,nSubj,nDays);
sdCursorAngle           = NaN(nBins,nSubj,nDays);
handAngle               = NaN(nBins,nSubj,nDays);
sdHandAngle             = NaN(nBins,nSubj,nDays);
reportAngle             = NaN(nBins,nSubj,nDays);
aimFixAngle             = NaN(nBins,nSubj,nDays);
implicitAngle           = NaN(nBins,nSubj,nDays);
implicitAngle_fromFix   = NaN(nBins,nSubj,nDays);
firstAimReportBin       = NaN(nSubj,1);
firstAimFixBin          = NaN(nSubj,1);

%% Loop over subjects

for s = 1 : nSubj
    
    load(subjFiles(s).name)
    disp(['Loaded ' subjFiles(s).name])
    subj{s} = ExpDetails(1).subjFolder;
    
    % get general experiment information
    if s==1
        iBin = Results.trialNo_bins(:,end);
        blockNo = Results.blockNo_bins(:,1);
        cursorRotation = Results.cursorRotation(iBin);
        rotBins = cursorRotation ~= 0;
        woBlock = find(ismember(lower(blockNames),'washout'));
        woBins = blockNo==woBlock;
    end
    
    %% Bin reaction time, cursor angle, hand angle, report angle, and fixation angle
    
    reportAngle_all = NaN(nBins,nDays);
    reportAngle_noOutliers = NaN(nBins,nDays);
    
    % loop over bins
    iBin = 1 : nTargets : nTrials+1;
    for b = 1 : length(iBin)-1
        % get all values in bin
        rt = Results.RT(iBin(b):iBin(b+1)-1,:);                 % reaction time
        ca = Results.hitAngle_cursor(iBin(b):iBin(b+1)-1,:);    % cursor hit angle
        ha = Results.hitAngle_hand(iBin(b):iBin(b+1)-1,:);      % hand hit angle
        ra_all = Results.explicitAngle(iBin(b):iBin(b+1)-1,:);  % reported aiming angle
        ra_no = Results.explicitAngle_outliersRemoved(iBin(b):iBin(b+1)-1,:);
        % get all aim fixation angles in bin
        fa = NaN;
        if rotBins(b)
            fa = Results.fixAngle_closestA(iBin(b):iBin(b+1)-1,:);
        elseif woBins(b)
            fa = Results.fixAngle_closestOppA(iBin(b):iBin(b+1)-1,:);
        end
        fa(:,sum(~isnan(fa))==1) = NaN;  % do not compute mean when there is only 1 value
        % calculate bin mean or median
        if meanOrMedian==1
            RT(b,s,:) = nanmean(rt);
            cursorAngle(b,s,:) = nanmean(ca);
            handAngle(b,s,:) = nanmean(ha);
            reportAngle_all(b,:) = nanmean(ra_all);
            reportAngle_noOutliers(b,:) = nanmean(ra_no);
            aimFixAngle(b,s,:) = nanmean(fa);
        elseif meanOrMedian==2
            RT(b,s,:) = nanmedian(rt);
            cursorAngle(b,s,:) = nanmedian(ca);
            handAngle(b,s,:) = nanmedian(ha);
            reportAngle_all(b,s,:) = nanmedian(ra_all);
            reportAngle_noOutliers(b,s,:) = nanmedian(ra_no);
            aimFixAngle(b,s,:) = nanmedian(fa);
        end
        % calculate bin standard deviation
        sdRT(b,s,:) = nanstd(rt);
        sdCursorAngle(b,s,:) = nanstd(ca);
        sdHandAngle(b,s,:) = nanstd(ha);
    end
    iBin = iBin(1:end-1)';
    
    % create separate aim fixation variables for rotation and washout?
    
    %% Check removed outliers in reported aiming angles
    
    keepOutliers = false;
    outlierDiff = ~isnan(reportAngle_noOutliers) & ...
        (reportAngle_all ~= reportAngle_noOutliers);
    if any(outlierDiff(:));
        outlier = ~isnan(Results.explicitAngle) & ...
            (Results.explicitAngle ~= Results.explicitAngle_outliersRemoved);
        [tr,day] = find(outlier);
        % plot reported values and removed outliers
        figure(fig1); clf
        plot([Results.explicitAngle_outliersRemoved],'o'); hold on
        for i = 1 : length(tr)
            plot(tr(i),Results.explicitAngle(tr(i),day(i)),'kx');
        end
        horline([-vmr,0,vmr])
        if nDays==2
            legend('Day1','Day2','outlier')
        end
        set(gca,'ytick',-180:45:180)
        xlabel('Trial'); ylabel('Angle (deg)')
        title('Reported aiming angles and removed outliers');
        keepOutliers = input('Discard(0) or keep(1) outliers? ');
    end
    if keepOutliers
        reportAngle(:,s,:) = reportAngle_all;
        keyboard
    else
        reportAngle(:,s,:) = reportAngle_noOutliers;
    end
    
    %% Calculate implicit angle
    
    implicitAngle(:,s,:) = handAngle(:,s,:)-reportAngle(:,s,:);
    implicitAngle_fromFix(:,s,:) = handAngle(:,s,:)-aimFixAngle(:,s,:);
    
    %% Trials and bins in which an aiming strategy was present
    
    % percentage of trials in which an aiming strategy was present
    reportTrial = ismember(Results.trialType,'report');
    for d = 1 : nDays
        for b = 1 : nBlocks
            % all trials
            currBlock = Results.blockNo(:,d)==b & Results.feedbackAndMTgood(:,d);
            strategy = Results.explicitAngle(currBlock,d)<targetZone(1) | ...
                Results.fixAngle_closestA(currBlock,d)<targetZone(1);
            percStrategy(b,s,d) = sum(strategy)/sum(currBlock)*100;
            strategy_opp = Results.explicitAngle(currBlock,d)>targetZone(2) | ...
                Results.fixAngle_closestOppA(currBlock,d)>targetZone(2);
            percStrategy_opp(b,s,d) = sum(strategy_opp)/sum(currBlock)*100;
            % report trials
            currBlock_report = Results.blockNo(:,d)==b & Results.feedbackAndMTgood(:,d) & reportTrial(:,d);
            if any(currBlock_report)
                reportStrategy = Results.explicitAngle(currBlock_report,d)<targetZone(1);
                percReportStrategy(b,s,d) = sum(reportStrategy)/sum(currBlock_report)*100;
                reportStrategy_opp = Results.explicitAngle(currBlock_report,d)>targetZone(2);
                percReportStrategy_opp(b,s,d) = sum(reportStrategy_opp)/sum(currBlock_report)*100;
            end
            % no-report trials, look at fixation angles
            currBlock_fix = Results.blockNo(:,d)==b & Results.feedbackAndMTgood(:,d) & ~reportTrial(:,d);
            if any(currBlock_fix)
                fixationStrategy = Results.fixAngle_closestA(currBlock_fix,d)<targetZone(1);
                percFixStrategy(b,s,d) = sum(fixationStrategy)/sum(currBlock_fix)*100;
                fixationStrategy_opp = Results.fixAngle_closestOppA(currBlock_fix,d)>targetZone(2);
                percFixStrategy_opp(b,s,d) = sum(fixationStrategy_opp)/sum(currBlock_fix)*100;
            end
        end
    end
    
    % first rotation bin in which an aiming strategy was present
    firstRotBin = find(rotBins,1);
    firstAimFixBin(s) = nanmax([find(aimFixAngle(firstRotBin:end,s,1)<targetZone(1),1) NaN]);
    firstAimReportBin(s) = nanmax([find(reportAngle(firstRotBin:end,s,1)<targetZone(1),1) NaN]);
    
    %% Plots
    
    if createPlots
        % TO DO: update code
        
        %% Plot binned angles - subplots
        
        %fig2 = scaledFigure(0.5+0.5*plotGaze+1*plotReport+0.5*plotRT,0.8*nDays);    % binned angles - subplots
        %fig2b = scaledFigure(1,0.8*nDays);                                          % binned angles - overlayed
        
        % dataToPlot = {Results.mnBinnedHitAngle_hand};
        % sdDataToPlot = {Results.sdBinnedHitAngle_hand};
        % figTitles = {'Learning'};
        % yLabels = {'Endpoint hand angle (deg)'};
        % legends = {[]};
        % if plotReport
        %     if ~isfield(Results,'mnBinnedExplicitAngle_outliersRemoved')
        %         Results.mnBinnedExplicitAngle_outliersRemoved = Results.binnedExplicitAngle;
        %     end
        %     Results.mnBinnedImplicitAngle_fix(~rotBins) = NaN;
        %     dataToPlot = [dataToPlot {Results.mnBinnedExplicitAngle_outliersRemoved,...
        %         cat(3,Results.mnBinnedImplicitAngle,Results.mnBinnedImplicitAngle_fix)}];
        %     sdDataToPlot = [sdDataToPlot {[],cat(3,[],[])}];
        %     figTitles = [figTitles {'Explicit learning','Implicit learning'}];
        %     yLabels = [yLabels {'Reported aim angle (deg)','Hand minus reported aim angle (deg)'}];
        %     legends = [legends {[],{'From report'}}];
        %     %legends = [legends {[],{'From report','From fixation'}}];
        % end
        % if plotGaze
        %     Results.mnBinnedFixAngle_closestA(~rotBins) = NaN;
        %     dataToPlot = [dataToPlot {cat(3,Results.mnBinnedFixAngle_closestA,...
        %         Results.mnBinnedFixAngle_closestOppA)}];
        %     sdDataToPlot = [sdDataToPlot {cat(3,Results.sdBinnedFixAngle_closestA,...
        %         Results.sdBinnedFixAngle_closestOppA)}];
        %     figTitles = [figTitles {'Preview fixation closest to hand target'}];
        %     yLabels = [yLabels {'Fixation angle (deg)'}];
        %     legends{end} = {'From report','From fixation'};
        %     legends = [legends {{['Closest to -' num2str(vmr)],['Closest to +' num2str(vmr)]}}];
        % end
        % if plotRT
        %     dataToPlot = [dataToPlot {Results.mnBinnedRT*1000}]; % in ms
        %     sdDataToPlot = [sdDataToPlot {Results.sdBinnedRT*1000}];
        %     figTitles = [figTitles {'Reaction time'}];
        %     yLabels = [yLabels {'Reaction time (ms)'}];
        %     legends = [legends {[]}];
        % end
        % colors = cat(3,colors,fadedColors); % concatenate for looping
        %
        % figure(fig2); clf;
        % nCol = length(dataToPlot);
        % for d = 1 : nDays
        %     for c = 1 : nCol
        %         % create subplot with correct size
        %         h(c) = subplot(nDays,nCol,d*nCol-nCol+c); hold on
        %         pb = pbaspect; pbaspect([pb(1) 0.95*pb(2) pb(3)]);
        %         h(c).Position(2) = h(c).Position(2)-0.025*pb(2);
        %         % set axes
        %         xlim([0 maxnTrials+1]); set(gca,'XTick',0:80:maxnTrials)
        %         if ~isempty(strfind(yLabels{c},'deg'));
        %             ylim([-95 50]); set(gca,'YTick',-90:45:45)
        %         elseif ~isempty(strfind(yLabels{c},'ms'));
        %             ylim([0 max(ExpDetails(1).timing.maxRT)*1000]);
        %         end
        %         yl = ylim;
        %         % color background when rotation is on
        %         for r = 1 : nRot
        %             x = iRotationOnOff(r*2-1:r*2,d);
        %             a = area([x; flipud(x)],[yl(1) yl(1) yl(2) yl(2)],'LineStyle','none');
        %             a.FaceColor = [0.95 0.95 0.95];
        %             plot(x,-[vmr vmr],'k') % draw line at hand target angle
        %         end
        %         plot([0 maxnTrials],[0 0],'k'); % draw line at 0
        %         % plot data
        %         for i = size(dataToPlot{c},3):-1:1
        %             p(i) = plot(iBin,dataToPlot{c}(:,d,i),'o','color',colors(d,:,i));
        %             if ~isempty(sdDataToPlot{c}(:,:,i))
        %                 errorb(iBin,dataToPlot{c}(:,d,i),sdDataToPlot{c}(:,d,i),'barwidth',0,'color',colors(d,:,i))
        %             end
        %         end
        %         hold off
        %         % add labels, title, and legend
        %         vertline(breakTrials(:,d),'g:')
        %         vertline([iNewBlock(:,d)-0.5; nTrials(d)+0.5],'k:')
        %         xlabel('Trial')
        %         ylabel(yLabels{c})
        %         if d==1
        %             title(figTitles{c})
        %             if ~isempty(legends{c})
        %                 legend([p(1) p(2)],legends{c},'location','southwest'); legend('boxoff')
        %             end
        %         end
        %     end
        % end
        % suplabel(['Binned angles - ' subjName],'t');
        %
        % % save
        % if savePlots
        %     saveFigAsPDF([saveToPath 'binnedAngles_' subjName],12)
        % end
        %
        % colors = colors(:,:,1); % reset
        %
        % %% Plot binned angles - single plot overlayed
        %
        % figure(fig2b); clf
        % for d = 1 : nDays
        %     % create subplot with correct size
        %     h(c) = subplot(nDays,1,d); hold on
        %     pb = pbaspect; pbaspect([pb(1) 0.95*pb(2) pb(3)]);
        %     h(c).Position(2) = h(c).Position(2)-0.025*pb(2);
        %     % set axes
        %     xlim([0 maxnTrials]); set(gca,'XTick',0:80:maxnTrials)
        %     ylim([-60 20]); yl = ylim; set(gca,'Ytick',-60:15:45)
        %     % color background when rotation is on
        %     for r = 1 : nRot
        %         x = iRotationOnOff(r*2-1:r*2,d);
        %         a = area([x; flipud(x)],[yl(1) yl(1) yl(2) yl(2)],'LineStyle','none');
        %         a.FaceColor = [0.95 0.95 0.95];
        %         plot(x,-[vmr vmr],'k') % draw line at hand target angle
        %     end
        %     plot([0 maxnTrials],[0 0],'k'); % draw line at 0
        %     % plot data
        %     pi = plot(iBin,Results.mnBinnedImplicitAngle(:,d),'o-','color',colors(5,:));
        %     ph = plot(iBin,Results.mnBinnedHitAngle_hand(:,d),'o-','color',colors(1,:));
        %     pe = plot(iBin,Results.mnBinnedExplicitAngle_outliersRemoved(:,d),'o-','color',colors(2,:));
        %     %pif = plot(iBin,Results.mnBinnedImplicitAngle_fix(:,d),'o-','color',fadedColors(5,:));
        %     pf = plot(iBin,Results.mnBinnedFixAngle_closestA(:,d),'o-','color',colors(4,:));
        %     hold off
        %     % add lables and title
        %     vertline(breakTrials(:,d),'g:')
        %     vertline([iNewBlock(:,d)-0.5; nTrials(d)+0.5],'k:')
        %     xlabel('Trial')
        %     ylabel('Angle (deg)')
        %     if d==1
        %         title(['Binned angles - ' subjName],'interpreter','none')
        %         legend([ph pe pi pf],{'Hand','Explicit','Implicit','Fixation'},'location','east');
        %         legend('boxoff')
        %     end
        % end
        %
        % if savePlots
        %     saveFigAsPDF([saveToPath 'binnedOverlayedAngles_' subjName],12)
        % end
        
    end
    
end % end of loop over subjects

%% Save data file

fileName = ['Results_' expName '_' datestr(now,'yyyymmdd') '.mat'];

% check if file does not exist yet
if exist([saveToPath fileName],'file') == 2
    disp(['A file named ' fileName ' already exists.'])
    overwrite = input('Do you want to overwrite it? Yes(1) or no(0): ');
else
    overwrite = 1;
end

% save file
if overwrite == 1
    save([saveToPath fileName],'subj','blockNames','blockNo','iBin','nTrials',...
        'cursorRotation','meanOrMedian','RT','handAngle','aimFixAngle',...
        'reportAngle','implicitAngle','implicitAngle_fromFix','firstAimReportBin','firstAimFixBin');
    disp(['Saved ' saveToPath fileName])
    
    % save copy of Matlab code used to process data
    mFilePath = mfilename('fullpath');
    saveCopyOfCode(mFilePath,saveToPath)
else
    disp('Results not been saved')
end

%% Plot all individual learning curves

iMidBin = iBin+3;
iRotationOnOff = find(diff(Results.cursorRotation(:,1))) + 0.5;
if mod(length(iRotationOnOff),2)
    iRotationOnOff = [iRotationOnOff; nTrials+0.5];
end
iNewBlock = find(diff(Results.blockNo(:,1))) + 0.5;

% colors
colors = brewermap(10,'Paired');
colors = colors([1,2,5,6,3,4,9,10],:); % blue, red, green, purple

% plot
scaledFigure(nDays*0.75,0.5*4);
subplotNo = 1;
for i = 1 : 4
    for d = 1 : nDays
        subplot(4,nDays,subplotNo); hold on
        xlim([0 iMidBin(end)+1]); set(gca,'XTick',0:80:nTrials)
        ylim([-2*vmr +vmr]); yl = ylim; set(gca,'Ytick',-180:45:180)
        % color background when rotation is on, works in R2015b but not in 2017a
        for r = 1 : length(iRotationOnOff)/2
            x = iRotationOnOff(r*2-1:r*2);
            a = area([x; flipud(x)],[yl(1) yl(1) yl(2) yl(2)],'LineStyle','none');
            a.FaceColor = [0.95 0.95 0.95];
        end
        plot([0 nTrials],[0 0],'k'); % draw line at zero
        plot([0 nTrials],[-vmr -vmr],'k'); % draw line at hand target angle
        % plot learning curves
        if i==1
            plot(iMidBin,handAngle(:,:,d),'.-','color',colors(1,:));
            plot(iMidBin,nanmean(handAngle(:,:,d),2),'.-','color',colors(2,:),'linewidth',2);
            figTitle = 'Hand angle';
        elseif i==2
            plot(iMidBin,reportAngle(:,:,d),'.-','color',colors(3,:));
            plot(iMidBin,nanmean(reportAngle(:,:,d),2),'.-','color',colors(4,:),'linewidth',2);
            figTitle = 'Report angle';
        elseif i==3
            if any(aimFixAngle(:))
                plot(iMidBin,aimFixAngle(:,:,d),'.-','color',colors(7,:));
                mAimFixAngle = nanmean(aimFixAngle(:,:,d),2);
                plot(iMidBin,mAimFixAngle,'.-','color',colors(8,:),'linewidth',2);
                nAimFixators = sum(~isnan(aimFixAngle(:,:,d)),2);
                str = ['max n = ' num2str(max(nAimFixators))];
                text(iMidBin(firstRotBin),0.5*vmr,str)
            end
            figTitle = 'Aimpoint fixation angle';
        elseif i==4
            % from report
            plot(iMidBin,implicitAngle(:,:,d),'.-','color',colors(5,:));
            mImplAngle = nanmean(implicitAngle(:,:,d),2);
            pr = plot(iMidBin,mImplAngle,'.-','color',colors(6,:),'linewidth',2);
            % from fixations
            implicitAngle_fromFix(~isnan(mImplAngle),:,d) = NaN;
            plot(iMidBin,implicitAngle_fromFix(:,:,d),'.--','color',colors(5,:));
            mImplAngle_fromFix = nanmean(implicitAngle_fromFix(:,:,d),2);
            pf = plot(iMidBin,mImplAngle_fromFix,'.--','color',colors(6,:),'linewidth',2);
            figTitle = 'Implicit angle';
            if d==1 && sum(~isnan(mImplAngle))>1 && sum(~isnan(mImplAngle_fromFix))>1
                legend([pr pf],{'from report','from fixation'},'Location','best','Orientation','horizontal');
            end
        end
        % axes and lables
        hold off
        vertline(breakTrials,':','color',[0.7 0.7 0.7])
        vertline(iNewBlock,'k:')
        if i==4
            xlabel('Trial')
        end
        if d==1
            ylabel('Angle (deg)')
        end
        if nDays==1
            title(figTitle)
        elseif nDays==2
            %title(['Day ' num2str(d) ' - ' figTitle])
            title([figTitle ' [day ' num2str(d) ']' ])
        end
        subplotNo = subplotNo+1;
    end
end
figTitle = [expName ' - group results (n=' num2str(nSubj) ')'];
suplabel(figTitle,'t');

% save figure
figName = [expName '_binnedAngles'];
if exist([saveToPath figName],'file') == 2 % check if file does not exist yet
    disp(['A figure named ' figName ' already exists.'])
    overwrite = input('Do you want to overwrite it? Yes(1) or no(0): ');
else
    overwrite = 1;
end
if overwrite == 1
    saveFigAsPDF(figName,12)
end
