function MakeFiguresForPaper(whichPanels, varargin)
% function MakeFiguresForPaper(whichPanels, ['dataDir' = 'originalData'])
%
% whichPanels = 1 => ROC/AUC distribution figure.
% whichPanels = 2 => Timecourse of binned AUC.

% In this script we make the ROC/AUC distribution figures.
p = inputParser;
p.addOptional('dataDir', 'originalData');
p.parse(varargin{:});

dataDir  = p.Results.dataDir;
currDir  = GetCurrentDirFromPathString(fileparts(mfilename('fullpath')));
figDir   = GetDataDirForFigure(6);

if (any(whichPanels==1))
  dataPath = fullfile(figDir, currDir, dataDir, 'Mfig_noduplicates.mat');
  data     = load(dataPath);
  M        = data.M;
  % Pick a few cell-odor pairs, along with the cell we show in the
  % raster. All indices past 81 need to have 1 subtracted from them, due
  % to the removal of the duplicate.
  dispPnIndsCmps = [46    2;
                    62    4;
                    82-1  5;
                    106-1 3;
                    122-1 6;
                    88-1  7]; % We show this cell in a raster

  pnRocInds = arrayfun(@(i) find(M.pnCells == dispPnIndsCmps(i,1) & M.pnCmps==dispPnIndsCmps(i,2)),1:size(dispPnIndsCmps,1));

  % We show all of the following cells in rasters
  dispKcIndsCmps = [84  5;
                    87  7;
                    29  6;
                    182 3;
                    183 5;
                    155 4];

  kcRocInds = arrayfun(@(i) find(M.kcCells == dispKcIndsCmps(i,1) & M.kcCmps==dispKcIndsCmps(i,2)),1:size(dispKcIndsCmps,1));

  % Get the true positives and false positives for the ROC curves
  pnFp = {M.RocPn(pnRocInds).x};
  pnTp = {M.RocPn(pnRocInds).y};

  kcFp = {M.RocKc(kcRocInds).x};
  kcTp = {M.RocKc(kcRocInds).y};

  % Start making the figures
  sfigure(FindFigureCreate('Figure 6E: PN and KC AUC')); clf; 
  set(gcf,'Color',[1 1 1],'NumberTitle','off');
  ResizeFigure(gcf,6,7,'inches');
  Q = ComputeSubplotPositions(2,2,[],0.1,0.015,0.1,0.00,0.05,0.1);

  % Roc curves
  subplotp(Q,1);
  h = PlotRocLikeKai(pnFp, pnTp);
  arrayfun(@(h) set(h,'Color',[0.5 0.75 1]), h);
  % Highlight the curve for one of the cells we've shown rasters for.
  set(h(end),'Color','b','LineWidth',2); 
  set(gca, 'xtick', [0 0.5 1.0], 'ytick', [0 0.5 1.0], 'FontSize', 12);
  set(gca, 'yticklabel', arrayfun(@num2str, get(gca, 'ytick'), 'UniformOutput', false));
  xlabel('FP rate', 'FontSize', 14);
  ylabel('TP rate', 'FontSize', 14);
  title('PNs', 'FontSize', 16);

  subplotp(Q,2);
  h = PlotRocLikeKai(kcFp, kcTp);
  arrayfun(@(h) set(h,'Color',[0.5 0.75 1]), h);
  % Highlight the curve for one of the cells we've shown rasters for.
  set(h(end),'Color','b','LineWidth',2);
  set(gca, 'xtick', [0 0.5 1.0], 'ytick', [0 0.5 1.0], 'yticklabel', [], 'FontSize', 12);
  xlabel('FP rate', 'FontSize', 14);
  title('KCs', 'FontSize', 16);

  % Auc distributions
  subplotp(Q,3);
  PlotAucDistributionLikeKai([M.RocPn.auc]);
  set(gca, 'xtick', [0 0.5 1.0], 'ytick', [0 25 50], 'FontSize', 12);
  xlabel('AUC',    'FontSize', 14);
  ylabel('#Cells', 'FontSize', 14);
  set(gca, 'yticklabel', arrayfun(@num2str, get(gca, 'ytick'), 'UniformOutput', false));

  subplotp(Q,4);
  PlotAucDistributionLikeKai([M.RocKc.auc]);
  set(gca, 'xtick', [0 0.5 1.0], 'ytick', [0 25 50], 'FontSize', 12);
  xlabel('AUC',    'FontSize', 14);

  arrayfun(@(i) set(subplotp(Q,i),'XColor',[0 0 0],'YColor',[0 0 0]),1:4);
end

if (any(whichPanels==2))
  % Plot the binned AUC figure.
  dataPath = fullfile(figDir, currDir, dataDir, 'Mbin50_noduplicates.mat');
  data     = load(dataPath);
  M        = data.M;

  pnAuc = reshape([M.RocPn.auc],size(M.RocPn));
  kcAuc = reshape([M.RocKc.auc],size(M.RocKc));
  sfigure(FindFigureCreate('Figure 6F: Binned PN and KC AUC')); clf;
  set(gcf,'Color',[1 1 1],'Resize','off','NumberTitle','off');
  ResizeFigure(gcf,14,7,'inches');
  hl = [];
  alpha = 0.01;
  [medPn, ciPn] = ComputePercentileAndConfidenceIntervalsFromData(pnAuc, alpha);
  [medKc, ciKc] = ComputePercentileAndConfidenceIntervalsFromData(kcAuc, alpha);

  binSize = 0.05;
  startTime = -1;
  t = (0:size(pnAuc,2)-1)*binSize + startTime;
  [hl(1), hp] = PlotMeanWithFilledCiBand(t, medPn, ciPn(2,:), ciPn(1,:), 'r', 2, 'r', 0.2);
  [hl(2), hp] = PlotMeanWithFilledCiBand(t, medKc, ciKc(2,:), ciKc(1,:), name2rgb('ForestGreen'), 2, name2rgb('ForestGreen'), 0.2);
  hl(3) = HorizontalLine(gca,0.5,'k','LineWidth',2);
  p1 = [];
  for i = 1:size(pnAuc,2)
    x = pnAuc(:,i); y = kcAuc(:,i);
    x = x(~isnan(x));
    y = y(~isnan(y));
    p1(i) = ranksum(x, y);
  end
  indSig = find(p1<0.01*2); % *2 to get a sided test.
  hl(5) = plot(t(indSig), 0.775*(indSig>0),'k*');
  ylim([0.4 0.8]);
  set(gca,'xtick',-1:3,'ytick',0.4:0.1:0.8);
 
  hl(4) = VerticalLine(gca,0,'k--','LineWidth',1);
  VerticalLine(gca,0.5,'k--','LineWidth',1);
  set(legend(hl,'PNs','KCs','Chance','Odor', 'p < 0.01','Location','NorthEast'),'box','off','Color','none');
  
  set(gca,'FontSize', 12);
  xlabel('Time (s)','FontSize', 14);
  ylabel('AUC','FontSize', 14);

  TightenAxesToFigure;
end