% function MakeFiguresForPaper(whichFigures)
%
% Makes the rand index figures for the paper.
function MakeFiguresForPaper(whichFigures, varargin)

p = inputParser;
p.addOptional('dataDir', 'originalData');
p.parse(varargin{:});

dataDir = p.Results.dataDir;
figDir  = GetDataDirForFigure(3);

thisDir     = GetCurrentDirFromPathString(fileparts(mfilename('fullpath')));
thisDirFull = fileparts(mfilename('fullpath'));

addpath(thisDirFull);

if (any(whichFigures==1))
  % Plot the the mean and sem of the ris for the different clustering methods.
  figId = 1;
  figName = 'RandIndicesComparingClusteringByDistancesWithClustByOdorAndClustByConc_MeanAndSem';

  dataFile = fullfile(figDir, thisDir, dataDir, 'Mfig.mat');
  data = load(dataFile);
  Mfig = data.Mfig{1};
  
  sfigure(FindFigureCreate('Figure 3H: Clustering Comparison')); clf; set(gcf,'Color',[1 1 1], 'NumberTitle', 'off');
  set(gcf,'Resize','off');
  ResizeFigure(gcf,8,8,'inches');
  b = linspace(0,1,31);
  
  uconcBl = mean(mean(Mfig.riByConcBl,1),3); % Mean across concs and kmeans
  uodorBl = mean(mean(Mfig.riByOdorBl,1),3);

  muConcBl = mean(uconcBl); % Mean across trials
  semConcBl = sem(uconcBl); % Sem across trials

  muOdorBl = mean(uodorBl);
  semOdorBl = sem(uodorBl);

  uconcResp = mean(mean(Mfig.riByConcResp,1),3); % Mean across concs and kmeans;
  uodorResp = mean(mean(Mfig.riByOdorResp,1),3); % Mean across concs and kmeans

  muConcResp = mean(uconcResp);
  semConcResp = sem(uconcResp);

  muOdorResp = mean(uodorResp);
  semOdorResp = sem(uodorResp);

  uchanceBl = mean(mean(mean(Mfig.riByChanceBl,1),3),4);  
  muChanceBl = mean(uchanceBl);
  seChanceBl = sem(uchanceBl);

  uchanceResp = mean(mean(mean(Mfig.riByChanceResp,1),3),4);  
  muChanceResp = mean(uchanceResp);
  seChanceResp = sem(uchanceResp);
  
  mu = [muConcBl muOdorBl muConcResp muOdorResp];
  lb = mu - [semConcBl semOdorBl semConcResp semOdorResp];
  ub = mu + [semConcBl semOdorBl semConcResp semOdorResp];
  colors = [1.00 0.75 0.75;
            0.75 0.75 1.00;
            1.00 0.00 0.00;
            0.00 0.00 1.00;];  
  hb = BarPlot(mu,lb,ub,0.4,colors,1); hold on;
  hch = plot([0.3 2.5], [muChanceBl muChanceBl],'--','Color',name2rgb('gray50'));
  plot([2.5 4.5], [muChanceResp muChanceResp],'--','Color',name2rgb('gray50'));
  xlim([0, 5]); 
  ylim([0.4 1]);
  set(gca,'xtick',[1.5 3.5],'ytick',[0.4:0.1:1],'xticklabel',{'Baseline', 'Response'},'FontSize',12);
  xlabel('Time Window','FontSize',14);
  ylabel('Rand Index', 'FontSize',14);
  title('Clustering of trajectories by odor vs. by concentration','FontSize',16);
  legend([hb(3:4) hch],'Clust. by dist. vs. by conc.','Clust. by dist. vs. by odor','Chance','Location','Northwest');
  disp('Mean and sem of the rand indices comparing clustering by distance with clustering by concentration (red), and for comparing clustering by distance with clustering by odor (blue), in the baseline [0.9 - 1.9] (pale), and response window ([2.1 3.1]), where the RI is first averagd over all 10 subsets of 3 of the 5 available concentrations.Ytick = [0.4:0.1:1]. Dashed lines are chance levels, computed both for baseline and response separately.');
  box on;
end
