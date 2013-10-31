function MakeFiguresForPaper(whichFigure, varargin)
% function MakeFiguresForPaper(whichFigure)
%
% Makes the correlation distance figures for the paper.
% whichFigure = 1 => Components and 1:1 mixture
% whichFigure = 2 => Summary
% whichFigure = 3 => Concentration series
% whichFigure = 4 => Summary
% whichFigure = 5 => Mixture morph
% whichFigure = 6 => Summary

p = inputParser;
p.addOptional('dataDir', 'originalData');
p.parse(varargin{:});

dataDir = p.Results.dataDir;
figDir  = GetDataDirForFigure(3);
currDir = GetCurrentDirFromPathString(fileparts(mfilename('fullpath')));
sourceDir = fullfile(figDir, currDir, dataDir);

switch whichFigure
 case 1 % Correlation distance for binary mixtures
  plotBinaryMixtureCorrelationDistance(fullfile(sourceDir, 'corrDistResults1to1Mixture.mat'), 'Figure 3B: Correlation distance');
 case 2
  plotBinaryMixtureCorrelationDistanceSummary(fullfile(sourceDir, 'corrDistResults1to1Mixture.mat'), 'Figure 3C: Summary operation');
 case 3 % Correlation distance for binary mixtures concentration series
  plotBinaryMixtureCorrelationDistance(fullfile(sourceDir, 'corrDistResultsConcentrationSeries.mat'),'Figure 3D: Concentration (correlation distances)');
 case 4
  plotBinaryMixtureCorrelationDistanceSummary(fullfile(sourceDir, 'corrDistResultsConcentrationSeries.mat'), 'Figure 3E: Concentration (summary)');
 case 5 % Correlation distance for binary mixtures
  plotBinaryMixtureCorrelationDistance(fullfile(sourceDir, 'corrDistResultsMixtureMorphs.mat'), 'Figure  3F: Odor morph (correlation distances)');
 case 6
  plotBinaryMixtureCorrelationDistanceSummary(fullfile(sourceDir, 'corrDistResultsMixtureMorphs.mat'), 'Figure  3F: Odor morph (summary)');
 otherwise
  error('Don''t know how to make figure %d.\n', whichFigure);
end

function plotBinaryMixtureCorrelationDistance(dataFile, figName)
data = load(dataFile);
sfigure(FindFigureCreate(figName)); ClearFigure;
set(gcf,'Resize','off','NumberTitle','off');
ResizeFigure(gcf,12,12,'inches');
ax = PlotCorrelationDistanceMatrices(data.D,0.005);

numPairs = size(data.whichPairs,1);
for i = 1:numPairs
  ttl = {sprintf('c%d',data.whichPairs(i,2)), sprintf('o%d', data.whichPairs(i,1))};
  title(ax(i), ttl,'FontSize',12);
end

for i = 1:numPairs
  ii  = (i-1)*numPairs+1;
  ttl = {sprintf('c%d',data.whichPairs(i,2)), sprintf('o%d', data.whichPairs(i,1))};
  ylabel(ax(ii), ttl,'FontSize',12);
end

function plotBinaryMixtureCorrelationDistanceSummary(dataFile, figName)
data = load(dataFile);
sfigure(FindFigureCreate(figName)); ClearFigure;
set(gcf,'Resize','off','NumberTitle','off');
ResizeFigure(gcf,12,12,'inches');
imagesc(data.Ds,[0 1]);

numPairs = size(data.whichPairs,1);
set(gca,'xtick',1:numPairs,'ytick',1:numPairs);

tickLabels = arrayfun(@(i) sprintf('c%d:o%d', data.whichPairs(i,2),data.whichPairs(i,1)), 1:numPairs, 'UniformOutput', false);
set(gca, 'xticklabel', tickLabels, 'ytickLabel', tickLabels,'xaxislocation','top','tickLength',[0 0]);

TightenAxesToFigure;
