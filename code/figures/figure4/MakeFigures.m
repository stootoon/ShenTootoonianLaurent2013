% function MakeFigures(['whichPanels'={'A','B','SA',...,'SG'}])
%
% Plots the panels in Figure 4 and S4. By default all panels will be
% plotted. Specific panels can be specified using the 'whichPanels'
% argument, provided with a cell array specifying the desired
% panels. Each element of the cell array should be the letter for the
% desired panel, prefixed by S if the panel is in the supplementary
% material.

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

