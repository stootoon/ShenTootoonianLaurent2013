% function MakeFigures(varargin)
% 
% Makes the panels for Figure 7 and S7.
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
% MakeFigures({'A','G'}, 'dataDir', 'recomputedData'); 

function MakeFigures(varargin)
whichFigure   = 7;
availablePanels = { 'A', 'G', 'H', 'J', 'K', 'SA'};

%% Parse the input arguments
[whichPanels, dataDir] = ParseArgsForMakeFigures(whichFigure, availablePanels, varargin);
if (isempty(whichPanels))
  return;
end

%% Setup function handles
MakeKcRasters      = FunctionHandleFromPath(fullfile(GetRootDir('figures'),'figure7','KcResponseLatency',                      'MakeFiguresForPaper.m'));
MakeKcRec          = FunctionHandleFromPath(fullfile(GetRootDir('figures'),'figure7','PnKcConnectivity2_noduplicates',         'MakeFiguresForPaper.m'));
MakePnRec          = FunctionHandleFromPath(fullfile(GetRootDir('figures'),'figure7','ComplexMixtureTrajectoryReconstruction', 'MakeFiguresForPaper.m'));
MakePnKcTimecourse = FunctionHandleFromPath(fullfile(GetRootDir('figures'),'figure7','PnKcResponseTimecourse_noduplicates',    'MakeFiguresForPaper.m'));

addpath([GetRootDir('figures'), '/figure7']);
addpath([GetRootDir('figures'), '/figure7/KcResponseLatency']);
addpath([GetRootDir('figures'), '/figure7/PnKcConnectivity2_noduplicates']);
addpath([GetRootDir('figures'), '/figure7/ComplexMixtureTrajectoryReconstruction']);
addpath([GetRootDir('figures'), '/figure7/PnKcResponseTimecourse_noduplicates']);

%% Make the plots
for i = 1:numel(whichPanels)
  panel = whichPanels{i};
  switch panel
   case 'A'
    MakeKcRasters(1);
   case {'B', 'C', 'D', 'E'}
    fprintf('Figure 7%s is not yet available.\n', panel);   
   case {'F', 'I'}
    fprintf('Figure 7%s is a schematic.\n', panel);
   case 'G'
    MakeKcRec(1, 'dataDir', dataDir);
   case 'H'
    MakeKcRec(2, 'dataDir', dataDir);
   case 'J'
    MakePnRec(1, 'dataDir', dataDir);
   case 'K'    
    MakePnRec(2, 'dataDir', dataDir);
   case 'SA'
    MakePnKcTimecourse();    
   otherwise
    error('Panel %d-%s does not exist.', whichFigure, panel);
  end
end

