function MakeFiguresForPaper(varargin)
% function MakeFiguresForPaper(['dataDir' = 'originalData'])

p = inputParser;
p.addOptional('dataDir','originalData');
p.addOptional('plotShuffles',false); % Can take a lot of RAM, so disable by default
p.parse(varargin{:});

plotShuffles = p.Results.plotShuffles;

dataDir = p.Results.dataDir;
currDir = GetCurrentDirFromPathString(fileparts(mfilename('fullpath')));
figDir  = GetDataDirForFigure(8);

dataFile = fullfile(figDir, currDir, dataDir, 'dataIdentityDecoding.mat');
PlotData = load(dataFile);

pnTimeBins = PlotData.Data.pnBinStarts;
kcTimeBins = PlotData.Data.kcBinStarts;

pnResults  = mean(PlotData.ResultsPn, 3); % Average over leave-one-out trials
kcResults  = mean(PlotData.ResultsKc, 3);

if (plotShuffles)
  % Load the Shuffled Pn Results
  dataShFile = fullfile(figDir, currDir, dataDir, 'dataIdentityDecodingPnRandomShuffles.mat');
  DataSh = load(dataShFile);
  numShuffles = numel(DataSh.Results);
  pnResultsSh = repmat(pnResults,[1 1 numShuffles]);
  for i = 1:numShuffles
    pnResultsSh(:,:,i) = mean(DataSh.Results{i}{1},3);
  end
end

t0 = 1.9;
t1 = 4.1;
whichPnTimeBins = find(pnTimeBins>=t0 & pnTimeBins<t1);
whichKcTimeBins = find(kcTimeBins>=t0 & kcTimeBins<t1);

odorOnsetTime = 2;
%% Identity Decoding
sfigure(FindFigureCreate('Figure 8A-C and S8A: Odor identification, categorization and generalization')); clf; 
set(gcf,'Resize','off','Color',[1 1 1],'NumberTitle','off');
ResizeFigure(gcf,9,9,'inches');
Q = ComputeSubplotPositions(3,2,[],0.08,0.01,0.02,0.05,0.075,0.1);

subplotp(Q,1);
tpn = pnTimeBins(whichPnTimeBins) - odorOnsetTime;
muPn = mean(pnResults(:,whichPnTimeBins));
semPn = std(pnResults(:,whichPnTimeBins))/sqrt(44);
hpn = PlotMeanWithFilledErrorBand(tpn, muPn, 3*semPn, 3*semPn, 'r', 2, 'r',0.2);

tkc = kcTimeBins(whichKcTimeBins) - odorOnsetTime;
muKc = mean(kcResults(:,whichKcTimeBins));
semKc = std(kcResults(:,whichKcTimeBins))/sqrt(44);
hkc = PlotMeanWithFilledErrorBand(tkc, muKc, 3*semKc, 3*semKc, name2rgb('ForestGreen'), 2, name2rgb('ForestGreen'),0.2);

if (plotShuffles)
  colSh = name2rgb('dodgerblue');
  % Plot the data for the shuffle with median peak performance
  muPnSh = squeeze(mean(pnResultsSh(:,whichPnTimeBins,:)));
  maxMuPnSh = max(muPnSh);
  whichShuffleToPlot = argmin(abs(maxMuPnSh - median(maxMuPnSh)));
  if (length(whichShuffleToPlot)>1) % If there's more than one closest to the median, pick the one with the highest max performance
    whichShuffleToPlot = whichShuffleToPlot(argmax(maxMuPnSh(whichShuffleToPlot),1));
  end

  muPnSh1 = muPnSh(:,whichShuffleToPlot);
  semPnSh1 = std(pnResultsSh(:,whichPnTimeBins,whichShuffleToPlot))/sqrt(44);  
  hsh = PlotMeanWithFilledErrorBand(tpn, muPnSh1, 3*semPnSh1, 3*semPnSh1,colSh,2,colSh,0.2);
end

ylim([0 1]);
xlim([t0 t1] - odorOnsetTime);
VerticalLine(gca,0.0, 'LineStyle', '--', 'Color', [0.5 0.5 0.5]);
VerticalLine(gca,0.5, 'LineStyle', '--', 'Color', [0.5 0.5 0.5]);
ylabel('Accuracy','FontSize', 14);
title('Identification','FontSize',16,'HorizontalAlignment','Right');

if (plotShuffles)
  set(legend([hpn hkc hsh], 'PNs','KCs','PN Subs.','Location','NorthEast'),'FontSize',10,'box','off','color','none');
else
  set(legend([hpn hkc], 'PNs','KCs','Location','NorthEast'),'FontSize',10,'box','off','color','none');
end

%% Plot the top n, median n, and bottom n performers
subplotp(Q,2);

hpn = plot(tpn, muPn, 'r','LineWidth',2); hold on;
hkc = plot(tkc, muKc, 'Color',name2rgb('ForestGreen'),'LineWidth',2); 

if (plotShuffles)
  numSh = 2;
  [foo, shOrd] = sort(maxMuPnSh);
  topSh = shOrd(end-numSh+1:end);
  botSh = shOrd(1:numSh);
  [foo, shOrdMedDist] = sort(abs(maxMuPnSh-median(maxMuPnSh)));
  medSh = shOrdMedDist(1:numSh);
  
  shCols = [name2rgb('royalblue');
            name2rgb('dodgerblue');
            name2rgb('deepskyblue1')];
  
  htop = plot(tpn,muPnSh(:,topSh),'Color',shCols(1,:),'LineWidth',2);
  hmed = plot(tpn,muPnSh(:,medSh),'Color',shCols(2,:),'LineWidth',2);
  hbot = plot(tpn,muPnSh(:,botSh),'Color',shCols(3,:),'LineWidth',2);
end  

ylim([0 1]);
xlim([t0 t1]-odorOnsetTime);

VerticalLine(gca,0.0, 'LineStyle', '--', 'Color', [0.5 0.5 0.5]);
VerticalLine(gca,0.5, 'LineStyle', '--', 'Color', [0.5 0.5 0.5]);
set(gca,'yticklabel', []);

if (plotShuffles)
  set(legend([htop(1) hmed(1) hbot(1)],'Best PN Subs.','Med. PN Subs.', 'Worst PN Subs.', 'location','northeast'),'FontSize',10,'color','none','box','off');
end

%% Category Decoding
dataShortFile = fullfile(figDir, currDir, dataDir, 'dataForClassificationShort.mat');
PlotData = load(dataShortFile);
t = PlotData.pnBinStarts - odorOnsetTime;

pnDataFile = fullfile(figDir, currDir, dataDir, 'dataCategoryDecodingPns.mat');
kcDataFile = fullfile(figDir, currDir, dataDir, 'dataCategoryDecodingKcs.mat');

PnData = load(pnDataFile);
KcData = load(kcDataFile);

pnResults = zeros(8, length(t));
kcResults = pnResults;

for i = 1:length('ABCDWXYZ')    
  pnResults(i,:) = mean(reshape(PnData.Results{i}{1},length(t),[]),2)';
  kcResults(i,:) = mean(reshape(KcData.Results{i}{1},length(t),[]),2)';
end

if (plotShuffles)
  pnShDataFile = fullfile(figDir, currDir, dataDir, 'dataCategoryDecodingPnRandomShuffles.mat');
  PnShData = load(pnShDataFile);
  
  numShuffles = numel(PnShData.Results)/8;
  pnShResults = zeros(8,length(t), numShuffles);
  
  for i = 1:numel(PnShData.Results)
    icmp = mod(i-1,8)+1;
    ish  = floor((i-1)/8)+1;
    pnShResults(icmp,:,ish) = mean(reshape(PnShData.Results{i}{1},length(t),[]),2)';
  end
end

% Now Plot
muPn = mean(pnResults); sdPn = std(pnResults); semPn = sdPn/sqrt(8);
muKc = mean(kcResults); sdKc = std(kcResults); semKc = sdKc/sqrt(8);

subplotp(Q,3); cla;
PlotMeanWithFilledErrorBand(t, muPn, 3*semPn, 3*semPn, 'r', 2, 'r',0.2);
PlotMeanWithFilledErrorBand(t, muKc, 3*semKc, 3*semKc, name2rgb('ForestGreen'), 2, name2rgb('ForestGreen'),0.2);
ylim([0.5 1]);
xlim([t0 t1]-odorOnsetTime);

if (plotShuffles)
  % Plot the shuffle closest to the median
  muPnSh = squeeze(mean(pnShResults));
  maxMuPnSh = max(muPnSh);
  whichShuffleToPlot = argmin(abs(maxMuPnSh - median(maxMuPnSh)));
  if (length(whichShuffleToPlot)>1)
    whichShuffleToPlot = whichShuffleToPlot(argmax(maxMuPnSh(whichShuffleToPlot),1));
  end
  muPnSh1 = mean(pnShResults(:,:,whichShuffleToPlot));
  sdPnSh1 = std(pnShResults(:,:,whichShuffleToPlot));
  semPnSh1= sdPnSh1/sqrt(8);
  
  PlotMeanWithFilledErrorBand(t, muPnSh1, 3*semPnSh1, 3*semPnSh1, colSh, 2, colSh,0.2);
end

VerticalLine(gca,0.0, 'LineStyle', '--', 'Color', [0.5 0.5 0.5]);
VerticalLine(gca,0.5, 'LineStyle', '--', 'Color', [0.5 0.5 0.5]);
ylabel('Accuracy','FontSize',14);
title('Categorization','FontSize',16,'HorizontalAlignment','Right');

%% In the next subplot, plot traces for the best, worst, and median shuffles
subplotp(Q,4);

plot(t, muPn, 'r','LineWidth',2); hold on;
plot(t, muKc, 'Color',name2rgb('ForestGreen'),'LineWidth',2); 

if (plotShuffles)
  [foo, shOrd] = sort(maxMuPnSh);
  topSh = shOrd(end-numSh+1:end);
  botSh = shOrd(1:numSh);
  [foo, shOrdMedDist] = sort(abs(maxMuPnSh-median(maxMuPnSh)));
  medSh = shOrdMedDist(1:numSh);
  
  plot(t, muPnSh(:,topSh),'Color',shCols(1,:),'LineWidth',2);
  plot(t, muPnSh(:,medSh),'Color',shCols(2,:),'LineWidth',2);
  plot(t, muPnSh(:,botSh),'Color',shCols(3,:),'LineWidth',2);
end

ylim([0.5 1]);
xlim([t0 t1]-odorOnsetTime);

VerticalLine(gca,0.0, 'LineStyle', '--', 'Color', [0.5 0.5 0.5]);
VerticalLine(gca,0.5, 'LineStyle', '--', 'Color', [0.5 0.5 0.5]);

set(gca,'yticklabel', []);

%% Generalization
PlotData = load(dataShortFile);
t = PlotData.pnBinStarts - odorOnsetTime;

pnDataFile = fullfile(figDir, currDir, dataDir, 'dataGeneralizationDecodingPns.mat');
kcDataFile = fullfile(figDir, currDir, dataDir, 'dataGeneralizationDecodingKcs.mat');
PnData = load(pnDataFile);
KcData = load(kcDataFile);

pnResults = zeros(8, length(t));
kcResults = pnResults;

for i = 1:length('ABCDWXYZ')    
  pnResults(i,:) = mean(reshape(nanmean(PnData.Results{i}{1}),length(t),[]),2)';
  kcResults(i,:) = mean(reshape(nanmean(KcData.Results{i}{1}),length(t),[]),2)';
end

if (plotShuffles)
  pnShDataFile = fullfile(figDir, currDir, dataDir, 'dataGeneralizationDecodingPnRandomShuffles.mat');
  PnShData     = load(pnShDataFile);
  numShuffles  = numel(PnShData.Results)/8;
  pnShResults  = zeros(8,length(t), numShuffles);
  
  for i = 1:numel(PnShData.Results)
    icmp = mod(i-1,8)+1;
    ish  = floor((i-1)/8)+1;
    pnShResults(icmp,:,ish) = mean(reshape(nanmean(PnShData.Results{i}{1}),length(t),[]),2)';
  end
end

% Now Plot
muPn = mean(pnResults); sdPn = std(pnResults); semPn = sdPn/sqrt(8);
muKc = mean(kcResults); sdKc = std(kcResults); semKc = sdKc/sqrt(8);

subplotp(Q,5); cla;
PlotMeanWithFilledErrorBand(t, muPn, 3*semPn, 3*semPn, 'r', 2, 'r',0.2);
PlotMeanWithFilledErrorBand(t, muKc, 3*semKc, 3*semKc, name2rgb('ForestGreen'), 2, name2rgb('ForestGreen'),0.2);
ylim([0.5 1]);
xlim([t0 t1]-odorOnsetTime);

if (plotShuffles)
  % Plot the shuffle closest to the median
  muPnSh = squeeze(mean(pnShResults));
  maxMuPnSh = max(muPnSh);
  whichShuffleToPlot = argmin(abs(maxMuPnSh - median(maxMuPnSh)));
  if (length(whichShuffleToPlot)>1)
    whichShuffleToPlot = whichShuffleToPlot(argmax(maxMuPnSh(whichShuffleToPlot),1));
  end
  muPnSh1 = mean(pnShResults(:,:,whichShuffleToPlot));
  sdPnSh1 = std(pnShResults(:,:,whichShuffleToPlot));
  semPnSh1= sdPnSh1/sqrt(8);
  
  PlotMeanWithFilledErrorBand(t, muPnSh1, 3*semPnSh1, 3*semPnSh1, colSh, 2, colSh,0.2);
end

VerticalLine(gca,0.0, 'LineStyle', '--', 'Color', [0.5 0.5 0.5]);
VerticalLine(gca,0.5, 'LineStyle', '--', 'Color', [0.5 0.5 0.5]);

ylabel('Accuracy','FontSize', 14);
xlabel('Time (s)','FontSize', 14);
title('Generalization','FontSize',16,'HorizontalAlignment','Right');

%% In the next subplot, plot traces for the best, worst, and median shuffles
subplotp(Q,6);

plot(t, muPn, 'r','LineWidth',2); hold on;
plot(t, muKc, 'Color',name2rgb('ForestGreen'),'LineWidth',2); 

if (plotShuffles)
  [foo, shOrd] = sort(maxMuPnSh);
  topSh = shOrd(end-numSh+1:end);
  botSh = shOrd(1:numSh);
  [foo, shOrdMedDist] = sort(abs(maxMuPnSh-median(maxMuPnSh)));
  medSh = shOrdMedDist(1:numSh);
  
  plot(t, muPnSh(:,topSh),'Color',shCols(1,:),'LineWidth',2);
  plot(t, muPnSh(:,medSh),'Color',shCols(2,:),'LineWidth',2);
  plot(t, muPnSh(:,botSh),'Color',shCols(3,:),'LineWidth',2);
end

ylim([0.5 1]);
xlim([t0 t1]-odorOnsetTime);
VerticalLine(gca,0.0, 'LineStyle', '--', 'Color', [0.5 0.5 0.5]);
VerticalLine(gca,0.5, 'LineStyle', '--', 'Color', [0.5 0.5 0.5]);
set(gca,'yticklabel', []);
