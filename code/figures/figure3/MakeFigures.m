% function MakeFigures(varargin)
% 
% Makes the panels for Figure 3 and S3.
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
whichFigure     = 3;
availablePanels = { 'A',  'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'SA','SB','SC','SD'};

%% Parse the input arguments
[whichPanels, dataDir] = ParseArgsForMakeFigures(whichFigure, availablePanels, varargin);
if (isempty(whichPanels))
  return;
end

%% Setup the function handles
thisDir               = fileparts(mfilename('fullpath'));

CorrDistPlotFunction  = FunctionHandleFromPath(fullfile(thisDir, 'CorrelationDistanceMatrices',        'MakeFiguresForPaper.m'));
PafPlotFunction       = FunctionHandleFromPath(fullfile(thisDir, 'BinaryMixtureTrajectoryProjections', 'MakeFiguresForPaper.m'));
RandIndexPlotFunction = FunctionHandleFromPath(fullfile(thisDir, 'BinaryMixtureTrajectoryClustering',  'MakeFiguresForPaper.m'));
MorphPlotFunction     = FunctionHandleFromPath(fullfile(thisDir, 'BinaryMixtureTrajectoryMorphs',      'MakeFiguresForPaper.m'));

%% Make the plots
for i = 1:numel(whichPanels)
  panel = whichPanels{i};
  switch panel
   case 'A'
    disp('Figure 3A: LLE plot not yet implemented.');
   case 'B'    
    CorrDistPlotFunction(1);
   case 'C'    
    CorrDistPlotFunction(2);
   case 'D'    
    CorrDistPlotFunction(2);
    disp('Figure 3D: LLE plot not yet implemented.');
   case 'E'    
    CorrDistPlotFunction(4);
    disp('Figure 3E: LLE plot not yet implemented.');
   case 'F'    
    CorrDistPlotFunction(6);    
    disp('Figure 3F: LLE plot not yet implemented.');
   case 'G'
    PafPlotFunction(1, 'dataDir', dataDir);
   case 'H'
    RandIndexPlotFunction(1, 'dataDir', dataDir);
   case 'I'
    MorphPlotFunction(1, 'dataDir', dataDir);
   case 'SA'
    PafPlotFunction(2,   'dataDir', dataDir);
   case 'SB'
    PafPlotFunction(3,   'dataDir', dataDir);
   case 'SC'
    MorphPlotFunction(2, 'dataDir', dataDir);
   case 'SD'    
    MorphPlotFunction(3, 'dataDir', dataDir);
    MorphPlotFunction(4, 'dataDir', dataDir);
   otherwise
    fprintf('Panel %d-%s not yet implemented.\n', whichFigure, panel);
  end
end

