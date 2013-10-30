% function ret = ValidatePanelSpecification(whichPanels, whichFigure, availablePanels)
%
% Helper function used to validate the panels specified for plotting. 
%
% INPUTS:
%
%   WHICHPANELS: A cell array specifying the desired panels using their
%   letters and with supplementary panels prefixed with 'S'. E.g. {'A','B',SA'};
%
%   WHICHFIGURE: An integer, the corresponding figure.
%
%   AVAILABLEPANELS: A cell array of available panels. 
%
% OUTPUT:
%
%   If inputs are validated, returns TRUE, otherwise generates an error.
function ret = ValidatePanelSpecification(whichPanels, whichFigure, availablePanels)
% Validate the inputs to this function
if (~iscell(whichPanels))
  error('Panels should be specified as a cell array.');
end

if (whichFigure<1 | whichFigure>8)
  error('whichFigure must be in [1-8].');
end

if (~iscell(availablePanels))
  error('availablePanels should be a cell array.');
end

% Validate the panel specification

ret = false;
for i = 1:numel(whichPanels)
  panel = whichPanels{i};
  validPanel = cellfun(@(x) isequal(upper(panel), upper(x)), availablePanels);
  if (all(~validPanel))
    fprintf('ERROR: Panel "%s" is not available for Figure %d.\n', panel, whichFigure);
    ListAvailablePanels(whichFigure, availablePanels);
    return;
  end
end
ret = true;