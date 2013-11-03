% function MakeFigures(varargin)
% 
% Makes the panels for Figure 4 and S4.
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
% MakeFigures({'A','B'}, 'dataDir', 'recomputedData'); 

function MakeFigures(varargin)
whichFigure = 4;
availablePanels = {'A', 'B', 'SA', 'SB', 'SC', 'SD', 'SE', 'SF', 'SG'};

%% Parse the input arguments
[whichPanels, dataDir] = ParseArgsForMakeFigures(whichFigure, availablePanels, varargin);
if (isempty(whichPanels))
  return;
end

%% Make the plots
thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);
for i = 1:numel(whichPanels)
  panel = whichPanels{i};
  switch panel
   case 'A'
    MakeFiguresForPaper(1);
   case 'B'
    MakeFiguresForPaper(2);
   case 'SA'
    MakeFiguresForPaper(2);
   case 'SB'
    MakeFiguresForPaper(2);
   case 'SC'
    MakeFiguresForPaper(2);
   case 'SD'    
    MakeFiguresForPaper(3);
   case 'SE'
    MakeFiguresForPaper(4);
   case 'SF'
    MakeFiguresForPaper(5);    
   case 'SG'
    MakeFiguresForPaper(6);
   otherwise
    error('Panel %d-%s do not exist.', whichFigure, panel);
  end
end

