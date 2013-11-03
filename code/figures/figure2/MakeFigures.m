% function MakeFigures(varargin)
% 
% Makes the panels for Figure 2 and S2.
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
whichFigure     = 2;
availablePanels = {'A', 'B', 'C', 'D', 'E', 'F', 'SA', 'SB', 'SC', 'SD', 'SE', 'SF'};

%% Parse the input arguments
[whichPanels, dataDir] = ParseArgsForMakeFigures(whichFigure, availablePanels, varargin);
if (isempty(whichPanels))
  return;
end
  
%% Setup the function handles
thisDir     = fileparts(mfilename('fullpath'));
MakeBinaryMixtureResponseLinearityFigures = FunctionHandleFromPath(fullfile(thisDir, 'BinaryMixtureResponseLinearity', 'MakeFiguresForPaper.m'));

%% Make the plots
plottedRasterExample = false;
for i = 1:numel(whichPanels)
  panel = whichPanels{i};
  switch panel
   case {'A','B','C','D'}
    if (~plottedRasterExample)
      MakeBinaryMixtureRasterExample;
      plottedRasterExample = true;
    end
   case 'E'
    MakeBinaryMixtureResponseLinearityFigures(4, 'dataDir', dataDir);
   case {'F', 'SA'}
    MakeBinaryMixtureResponseLinearityFigures(1, 'dataDir', dataDir);
   case 'SB'
    MakeBinaryMixtureResponseLinearityFigures(2, 'dataDir', dataDir);
   case 'SC'
    MakeBinaryMixtureResponseLinearityFigures(3, 'dataDir', dataDir);
   case 'SD'    
    MakeBinaryMixtureResponseLinearityFigures(5, 'dataDir', dataDir);
   case 'SE'
    MakeBinaryMixtureResponseLinearityFigures(6, 'dataDir', dataDir);
   case 'SF'       
    MakeBinaryMixtureResponseLinearityFigures(7, 'dataDir', dataDir);
   otherwise
    fprintf('Panel %d-%s do not exist.', whichFigure, panel);
  end
end
