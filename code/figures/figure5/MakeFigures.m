% function MakeFigures(['whichPanels' = {'D','E','F','SA',...,'SD'}], ['dataDir' = 'originalData'])
%
% Plots the panels in Figure 5 and S5. By default all panels will be
% plotted. Specific panels can be specified using the 'whichPanels'
% argument, provided with a cell array specifying the desired
% panels. Each element of the cell array should be the letter for the
% desired panel, prefixed by S if the panel is in the supplementary
% material.

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

