% function MakeFigures(varargin)
% 
% Makes the panels for Figure 8 and S8.
%
% Usage:
%
% List available panels:        MakeFigures();
% Plot one panel e.g. A:        MakeFigures('A');
% Plot multiple panels:         MakeFigures({'A','SA'});
% Plot all panels:              MakeFigures('all');
%
% If you've run ProcessData to recompute the data for the
% figures, you can use its results by specifying 'dataDir':
%
% Plot multiple panels, using the recomputed data: 
% MakeFigures({'A','C'}, 'dataDir', 'recomputedData'); 

function MakeFigures(varargin)
whichFigure   = 8;
availablePanels = {'A', 'B', 'C', 'SA', 'SB', 'SC'};

%% Parse the input arguments
[whichPanels, dataDir] = ParseArgsForMakeFigures(whichFigure, availablePanels, varargin);
if (isempty(whichPanels))
  return;
end

%% Setup function handles
MakePerComp = FunctionHandleFromPath(fullfile(GetRootDir('figures'),'figure8','generalization_noduplicates', 'MakeFiguresForPaper.m'));
MakeOverall = FunctionHandleFromPath(fullfile(GetRootDir('figures'),'figure8','classification_noduplicates', 'MakeFiguresForPaper.m'));

addpath([GetRootDir('figures'), '/figure8']);
addpath([GetRootDir('figures'), '/figure8/generalization_noduplicates']);
addpath([GetRootDir('figures'), '/figure8/classification_noduplicates']);

%% Make the plots
% The shuffles files are big and loading them takes a lot of time and
% RAM. Set this to true if you want to plot them anyway!
plotShuffles = false; 

plottedPerformance = false;
for i = 1:numel(whichPanels)
  panel = whichPanels{i};
  switch panel
   case {'A', 'B', 'C', 'SA'}
    if (~plottedPerformance)
      MakeOverall('dataDir', dataDir, 'plotShuffles', plotShuffles);
      plottedPerformance = true;
    end
   case {'D', 'E','F','G','H','I'}
    fprintf('Panel 8%s is a schematic.\n', panel);
   case 'SB'
    MakePerComp(1, 'dataDir', dataDir);
   case 'SC'
    MakePerComp(2, 'dataDir', dataDir);
   otherwise
    error('Panel %d-%s does not exist.', whichFigure, panel);
  end
end

