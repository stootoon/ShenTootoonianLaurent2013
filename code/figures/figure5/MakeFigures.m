% function MakeFigures(varargin)
% 
% Makes the panels for Figure 5 and S5.
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
whichFigure     = 5;
availablePanels = {'A', 'B', 'C', 'D', 'E', 'F', 'SA', 'SB', 'SC', 'SD'};

%% Parse the input arguments
[whichPanels, dataDir] = ParseArgsForMakeFigures(whichFigure, availablePanels, varargin);
if (isempty(whichPanels))
  return;
end
  
%% Setup the function handles
thisDir = fileparts(mfilename('fullpath')); 
addpath(thisDir);

MakeComplexMixtureGlobalTrajectoryDistance = FunctionHandleFromPath(fullfile(GetRootDir('figures'),'figure5','ComplexMixtureGlobalTrajectoryDistanceWithOverlap_noduplicates','MakeFiguresForPaper.m'));
MakeComplexMixtureTransientProximity       = FunctionHandleFromPath(fullfile(GetRootDir('figures'),'figure5','ComplexMixtureTransientProximity_noduplicates',                 'MakeFiguresForPaper.m'));
MakeComplexMixturePerBinTrajectoryDistance = FunctionHandleFromPath(fullfile(GetRootDir('figures'),'figure5','ComplexMixturePerBinTrajectoryDistanceWithOverlap_noduplicates','MakeFiguresForPaper.m'));
MakeCorrelationDistanceSummary             = FunctionHandleFromPath(fullfile(GetRootDir('figures'),'figure5','CorrelationDistanceMatrices','MakeFiguresForPaper.m'));

%% Now make the plots
for i = 1:numel(whichPanels)
  panel = whichPanels{i};
  switch panel
   case 'A'
    MakeCorrelationDistanceSummary(2);
    disp('Figure 5A: LLE not yet implemented.');
   case 'B'
    MakeCorrelationDistanceSummary(4);
    disp('Figure 5B: LLE not yet implemented.');
   case 'C'
    MakeCorrelationDistanceSummary(6);
    disp('Figure 5C: LLE not yet implemented.');
   case 'D'
    MakeComplexMixtureGlobalTrajectoryDistance(1, 'dataDir', dataDir);
   case 'E'
    MakeComplexMixtureTransientProximity(1, 'dataDir', dataDir);
   case 'F'    
    MakeComplexMixturePerBinTrajectoryDistance(1, 'dataDir', dataDir);
   case 'SA'
    MakeComplexMixtureGlobalTrajectoryDistance(2, 'dataDir', dataDir);
   case 'SB'
    MakeComplexMixtureGlobalTrajectoryDistance(3, 'dataDir', dataDir);
   case 'SC'
    MakeComplexMixturePerBinTrajectoryDistance(2, 'dataDir', dataDir);
   case 'SD'
    MakeComplexMixturePerBinTrajectoryDistance(3, 'dataDir', dataDir);
   otherwise
    error('Panel %d-%s do not exist.', whichFigure, panel);
  end
end

