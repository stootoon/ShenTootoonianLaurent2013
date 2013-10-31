% function MakeFigures(['whichPanels'={'A',...,'I','SA',...,'SC'}])
%
% Plots the panels in Figure 8 and S8. By default all panels will be
% plotted. Specific panels can be specified using the 'whichPanels'
% argument, provided with a cell array specifying the desired
% panels. Each element of the cell array should be the letter for the
% desired panel, prefixed by S if the panel is in the supplementary
% material.

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

