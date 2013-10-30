% function MakeFigures(varargin)
% 
% 
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
