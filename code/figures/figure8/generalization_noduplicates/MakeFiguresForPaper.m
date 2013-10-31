function MakeFiguresForPaper(whichPanels, varargin)
% MakeFiguresForPaper(whichPanels, ['dataDir' = 'originalData'])
%
% whichPanels = 1 => S8B: Time course of categorization performance per component.
% whichPanels = 2 => S8C: Time course of generalization performance per component.

p = inputParser;
p.addOptional('dataDir', 'originalData');
p.parse(varargin{:});

dataDir = p.Results.dataDir;
currDir = GetCurrentDirFromPathString(fileparts(mfilename('fullpath')));
figDir  = GetDataDirForFigure(8);

LOADF = @(fileName) load(fullfile(figDir, currDir, dataDir, fileName));

cmps = 'ABCDWXYZ';
if (any(whichPanels==1))
  data = LOADF('CatResults.mat'); % Categorization results, computed Kai's way
  CatResults = data.CatResults;
  numBins = size(CatResults{1}{2},1);
  pnMean = nan(numBins,8);
  kcMean = nan(numBins,8);
  pnStd  = pnMean;
  pnSe   = pnStd;
  
  kcStd  = kcMean;
  kcSe   = kcStd;
  
  t = (0:numBins-1)*0.025;
  maxBins = 77; % To match Kai's plots

  sfigure(FindFigureCreate('Figure S8B: Time course of categorization performance per odor component'));
  clf; set(gcf,'Color',[1 1 1],'NumberTitle','off'); ResizeFigure(gcf,12,6,'inches');
  Q = ComputeSubplotPositions(2,4,[],0.05,0.01,0.02,0.01,0.1,0.05);
  for i = 1:8
    [pnMean(:,i), pnStd(:,i), pnSe(:,i)] = computeCategorizationPerformance(CatResults{i}{2});
    [kcMean(:,i), kcStd(:,i), kcSe(:,i)] = computeCategorizationPerformance(CatResults{i}{3}); 
    subplotp(Q,i);
    plotSiFigure(t(1:maxBins), pnMean(1:maxBins,i), 3*pnSe(1:maxBins,i), kcMean(1:maxBins,i), 3*kcSe(1:maxBins,i));
    text(1.8,0.95, cmps(i),'FontSize',14,'HorizontalAlignment','left','VerticalAlignment','bottom');
    if (i==5)
      set(gca,'xticklabel',arrayfun(@num2str,get(gca,'xtick'),'UniformOutput', false));
      set(gca,'yticklabel',arrayfun(@num2str,get(gca,'ytick'),'UniformOutput', false));
      set(gca,'FontSize', 12);
      xlabel('Time (s)','FontSize',14);
      ylabel('Accuracy','FontSize',14);
    end
  end

end

if (any(whichPanels==2))
  data = LOADF('GenResults.mat'); % Generalization results, computed Kai's way
  GenResults = data.GenResults;
  numBins = size(GenResults{1}{2},2);
  pnMean = nan(numBins,8);
  kcMean = nan(numBins,8);
  pnStd  = pnMean;
  pnSe   = pnStd;
  
  kcStd  = kcMean;
  kcSe   = kcStd;
  
  t = (0:numBins-1)*0.025;
  maxBins = 77; % To match Kai's plots

  sfigure(FindFigureCreate('Figure S8C: Time course of generalization performance per odor component'));
  clf; set(gcf,'Color',[1 1 1],'NumberTitle','off'); ResizeFigure(gcf,12,6,'inches');

  Q = ComputeSubplotPositions(2,4,[],0.05,0.01,0.02,0.01,0.1,0.05);
  for i = 1:8
    [pnMean(:,i), pnStd(:,i), pnSe(:,i)] = computeGeneralizationPerformance(GenResults{i}{2});
    [kcMean(:,i), kcStd(:,i), kcSe(:,i)] = computeGeneralizationPerformance(GenResults{i}{3}); 
    subplotp(Q,i);
    plotSiFigure(t(1:maxBins), pnMean(1:maxBins,i), 3*pnSe(1:maxBins,i), kcMean(1:maxBins,i), 3*kcSe(1:maxBins,i));
    text(1.8,0.95, cmps(i),'FontSize',14,'HorizontalAlignment','left','VerticalAlignment','bottom');
    if (i==5)
      set(gca,'xticklabel',arrayfun(@num2str,get(gca,'xtick'),'UniformOutput', false));
      set(gca,'yticklabel',arrayfun(@num2str,get(gca,'ytick'),'UniformOutput', false));
      set(gca,'FontSize', 12);
      xlabel('Time (s)','FontSize',14);
      ylabel('Accuracy','FontSize',14);
    end
  end
end

function plotSiFigure(t,pnMean, pnSem, kcMean, kcSem)

PlotMeanWithFilledErrorBand(t, pnMean, pnSem,pnSem,'r',2,'r',0.25);
PlotMeanWithFilledErrorBand(t, kcMean, kcSem,kcSem,name2rgb('ForestGreen'),2,name2rgb('ForestGreen'),0.25);  

maxPn = max(pnMean); tmaxPn = t(argmax(pnMean,1));
maxKc = max(kcMean); tmaxKc = t(argmax(kcMean,1));
plot(tmaxPn, maxPn,'ro','MarkerSize',10,'MarkerFaceColor','r');
plot(tmaxKc, maxKc,'o','Color', name2rgb('ForestGreen'),'MarkerSize',10,'MarkerFaceColor',name2rgb('ForestGreen'));

ylim([0.5 1]);
xlim([0 2]);

VerticalLine(gca,0.1, 'b');
VerticalLine(gca,0.6, 'b');

set(gca,'xtick',[0 1 2]);
set(gca,'ytick',[0:0.1:1]);
set(gca,'xticklabel',[],'yticklabel',[]);
set(gca,'linewidth',0.1);
box on;
axis square;

function [rm, rsd, rse] = computeGeneralizationPerformance(R)
R = squeeze(nanmean(R)); % Average over the odors left out
R = reshape(R,size(R,1),[]); % bins x [bootstraps x valences]
rm = mean(R,2); % Average over bootstraps and valences
rsd = std(R,[],2); 
rse = rsd/sqrt(size(R,2));


function [rm, rsd, rse] = computeCategorizationPerformance(R)
R = permute(R,[3 1 2 4]); % Permute to match the format of the Generalization results
R = squeeze(nanmean(R)); % Average over the odors left out
R = reshape(R,size(R,1),[]); % bins x [bootstraps x valences]
rm = mean(R,2); % Average over bootstraps and valences
rsd = std(R,[],2); 
rse = rsd/sqrt(size(R,2));
