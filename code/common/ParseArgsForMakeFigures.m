function [whichPanels, dataDir] = ParseArgsForMakeFigures(whichFigure, availablePanels, vargs)
% [whichPanels, dataDir] = ParseArgsForMakeFigures(whichFigure, availablePanels, vargs)
%
% Parses the arguments to a MAKEFIGURES call, validates the panel
% selection, and returns the list of panels selected and the desired
% data directory.

% Parse the arguments

whichPanels = '';
dataDir     = 'originalData';
nargs = numel(vargs);
switch(nargs)
 case 0 % No input args, list the available panels.
  ListAvailablePanels(whichFigure, availablePanels);
  return;
 case 1 % One input argument, take it to the panels to be plotted
  arg = vargs{1};
  if (~iscell(arg))
    if (isequal(lower(arg),'all'))
      whichPanels = availablePanels;
    else
      whichPanels = {vargs{1}};
    end
  else
    whichPanels = vargs{1};
  end
 otherwise
  p = inputParser;
  p.addOptional('whichPanels', availablePanels)
  p.addOptional('dataDir',     'originalData');
  p.parse(vargs{:});
  whichPanels = p.Results.whichPanels;
  dataDir     = p.Results.dataDir;
end

% Validate the panel specfiication
validSpecs = ValidatePanelSpecification(whichPanels, whichFigure, availablePanels);
if (~validSpecs)
  error('Invalid panel specification.');
end

