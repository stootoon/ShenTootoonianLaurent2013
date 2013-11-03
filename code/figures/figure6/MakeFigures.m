% function MakeFigures(varargin)
% 
% Makes the panels for Figure 6 and S6.
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
whichFigure   = 6;
availablePanels = {'A', 'B', 'C', 'D', 'E', 'F', 'SA', 'SB','SC','SD'};

%% Parse the input arguments
[whichPanels, dataDir] = ParseArgsForMakeFigures(whichFigure, availablePanels, varargin);
if (isempty(whichPanels))
  return;
end

%% Setup the function handles
MakePnKcAuc      = FunctionHandleFromPath(fullfile(GetRootDir('figures'),'figure6','PnKcAuc_noduplicates','MakeFiguresForPaper.m'));
MakeResponsivity = FunctionHandleFromPath(fullfile(GetRootDir('figures'),'figure6','PnAndKcPromiscuity_noduplicates','MakeFiguresForPaper.m'));

addpath([GetRootDir('figures'),'/figure6']);
addpath([GetRootDir('figures'),'/figure6/PnKcAuc_noduplicates']);
addpath([GetRootDir('figures'),'/figure6/PnAndKcPromiscuity_noduplicates']);

%% Now make the plots
for i = 1:numel(whichPanels)
  panel = whichPanels{i};
  switch panel
   case 'A'
    sfigure(FindFigureCreate('Figure 6A: PN 1 Odor Responses')); clf;
    set(gcf,'NumberTitle','off');
    pnSpt = LoadTocSpikeTimes('rawpn');
    PlotOdorResponseRastersForCell(pnSpt, 87, [1.7 4.2], 'F');
    ResizeFigure(gcf,12,8,'inches');
    set(gcf,'Name', 'Figure 6A: PN 1 Odor Responses', 'NumberTitle', 'off');
   
   case 'B'
    sfigure(FindFigureCreate('Figure 6B: KC 1 Odor Responses')); clf;
    set(gcf, 'NumberTitle', 'off');
    kcSpt = LoadTocSpikeTimes('rawkc');
    PlotOdorResponseRastersForCell(kcSpt, 155, [1.7 4.2], 'D', 'spikeColor','k');
    ResizeFigure(gcf,12,8,'inches');
    set(gcf,'Name', 'Figure 6B: KC 1 Odor Responses', 'NumberTitle', 'off');    
   
   case 'C'
    kcSpt = LoadTocSpikeTimes('rawkc');
    
    sfigure(FindFigureCreate('Figure 6C: KC 2 Odor Responses')); clf;
    set(gcf, 'NumberTitle', 'off');
    PlotOdorResponseRastersForCell(kcSpt, 84, [1.9 3.2], 'W', 'spikeColor','k');
    ResizeFigure(gcf,12,8,'inches');

    sfigure(FindFigureCreate('Figure 6C: KC 3 Odor Responses')); clf;
    set(gcf, 'NumberTitle', 'off');    
    PlotOdorResponseRastersForCell(kcSpt, 87, [1.9 3.2], 'Y', 'spikeColor','k');
    ResizeFigure(gcf,12,8,'inches');
    
   case 'D'
    kcSpt = LoadTocSpikeTimes('rawkc');
    
    sfigure(FindFigureCreate('Figure 6D: KC 4 Odor Responses')); clf;
    set(gcf, 'NumberTitle', 'off');
    PlotOdorResponseRastersForCell(kcSpt, 29, [1.9 3.9], 'X', 'spikeColor','k');
    ResizeFigure(gcf,12,8,'inches');

    sfigure(FindFigureCreate('Figure 6D: KC 5, 6 Odor Responses')); clf;
    set(gcf, 'NumberTitle', 'off');
    PlotOdorResponseRastersForCell(kcSpt, 183, [1.95 3.5], 'F', 'spikeColor','r','spikeWidth',1);
    PlotOdorResponseRastersForCell(kcSpt, 182, [1.95 3.5], 'F', 'spikeColor','b','spikeWidth',1,'clearFigure',false);
    ResizeFigure(gcf,12,8,'inches');    
   case 'E'    
    MakePnKcAuc(1, 'dataDir', dataDir);
   case 'F'
    MakePnKcAuc(2, 'dataDir', dataDir);
   case 'SA'
    MakeResponsivity(3); % Responsivities are computed directly from the raw data rather than processed data, so don't need 'dataDir'.
   case 'SB'
    MakeResponsivity(4);
   case 'SC'
    MakeResponsivity(1);
   case 'SD'
    MakeResponsivity(2);
   otherwise
    error('Panel %d-%s do not exist.', whichFigure, panel);
  end
end

