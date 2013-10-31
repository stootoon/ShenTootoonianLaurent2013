function MakeFiguresForPaper(whichFigure, varargin)
% function MakeFiguresForPaper(whichFigure, varargin)
% 
% whichFigure=1,2  => Figure 5A: Single odors and summary
% whichFigure=3,4  => Figure 5B: W-family and summary
% whichFigure=5,6  => Figure 5C: Two families and summary
% whichFigure=7,8  => Figure 5-: Partial overlaps
% whichFigure=9,10 => Figure 5E: Transient overlap.

p = inputParser;
p.addOptional('dataDir', 'originalData');
p.parse(varargin{:});

dataDir = p.Results.dataDir;
figDir  = GetDataDirForFigure(5);
thisDir = GetCurrentDirFromPathString(fileparts(mfilename('fullpath')));

sourceDir = fullfile(figDir, thisDir, dataDir);

switch whichFigure
 case 1 % Correlation distance for single components.
  plotBinaryMixtureCorrelationDistance(fullfile(sourceDir, 'corrDistResultsSingleComponents.mat'), 'Figure 5A: Single odors');
 case 2 % Summary
  plotBinaryMixtureCorrelationDistanceSummary(fullfile(sourceDir, 'corrDistResultsSingleComponents.mat'), 'Figure 5A: Single odors (summary)');
 case 3 % Correlation distance for W series
  plotBinaryMixtureCorrelationDistance(fullfile(sourceDir, 'corrDistResultsWseries.mat'), 'Figure 5B: W-family');
 case 4 % Summary 
  plotBinaryMixtureCorrelationDistanceSummary(fullfile(sourceDir, 'corrDistResultsWseries.mat'), 'Figure 5B: W-family (summary)');
 case 5 % Correlation distance for A vs W series
  plotBinaryMixtureCorrelationDistance(fullfile(sourceDir, 'corrDistResultsWfamilyAfamily.mat'), 'Figure 5C: Two families');
 case 6 % Summary
  plotBinaryMixtureCorrelationDistanceSummary(fullfile(sourceDir, 'corrDistResultsWfamilyAfamily.mat'), 'Figure 5C: Two families (summary)');
 case 7 % Correlation distance for partial overlaps
  plotBinaryMixtureCorrelationDistance(fullfile(sourceDir, 'corrDistResultsPartialOverlaps.mat'), 'Figure 5-: Partial overlaps');
 case 8 % Summary
  plotBinaryMixtureCorrelationDistanceSummary(fullfile(sourceDir, 'corrDistResultsPartialOverlaps.mat'), 'Figure 5-: Partial overlaps (summary)');
 case 9 % Correlation distance for transient overlaps
  plotBinaryMixtureCorrelationDistance(fullfile(sourceDir, 'corrDistResultsTransientOverlaps.mat'), 'Figure 5E: Transient overlap ');
 case 10 % Summary
  plotBinaryMixtureCorrelationDistanceSummary(fullfile(sourceDir, 'corrDistResultsTransientOverlaps.mat'), 'Figure 5E: Transient overlap (summary)');
 otherwise
  error('Don''t know how to make figure %d.\n', whichFigure);
end

function plotBinaryMixtureCorrelationDistance(dataFile, figName)
data = load(dataFile);
sfigure(FindFigureCreate(figName)); ClearFigure;
set(gcf,'Resize','off','NumberTitle','off');
ResizeFigure(gcf,12,12,'inches');
ax = PlotCorrelationDistanceMatrices(data.D,0.005);

numOdors = numel(data.whichOdors);
for i = 1:numOdors
  title(ax(i), data.whichOdors{i},'FontSize',12);
end

for i = 1:numOdors
  ii  = (i-1)*numOdors+1;
  ylabel(ax(ii), data.whichOdors{i},'FontSize',12);
end

function plotBinaryMixtureCorrelationDistanceSummary(dataFile, figName)
data = load(dataFile);
sfigure(FindFigureCreate(figName)); ClearFigure;
set(gcf,'Resize','off','NumberTitle','off');
ResizeFigure(gcf,12,12,'inches');
imagesc(data.Ds,[0 1]);

numOdors = numel(data.whichOdors);
set(gca, 'xtick', 1:numOdors, 'ytick', 1:numOdors);
set(gca, 'xticklabel', data.whichOdors, 'ytickLabel', data.whichOdors, 'xaxislocation', 'top', 'tickLength', [0 0]);

TightenAxesToFigure;
