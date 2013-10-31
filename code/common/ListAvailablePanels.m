function ListAvailablePanels(whichFigure, panelsList)
% ListAvailablePanels(whichFigure, panelsList)
%
% Lists the available panels specified, and provides some usage information.
fprintf('Panels available for Figure %d: ', whichFigure);
cellfun(@(s) fprintf('%s ', s), panelsList);
fprintf('\n\n');
fprintf('USAGE:\n');
fprintf('Plot one panel e.g. ''A'':\t\tMakeFigures(''A''); \n');
fprintf('Plot several panels e.g. ''A'',''SB'':\tMakeFigures({''A'',''SB''}); \n');
fprintf('Plot all panels:\t\t\tMakeFigures(''all''); \n');
