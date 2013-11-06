function SpikeTimesBrowser(varargin)
% SpikeTimesBrowser(varargin)
%
% Plots the spike times of single cells in the binary mixtures or
% complex mixtures experiments as rasters. Binary mixtures rasters are
% plotted as in Figure 2A-D, and complex mixtures rasters are plotted
% as in Figure 6A. Once the figure has focus (by clicking on it), the
% left and right arrow keys can be used to progress through the cells
% in the dataset.
%
% SpikeTimesBrowser() called with now arguments plots the PNs in the
% binary mixture experiments. Other cells and settings can be
% specified using the following options, supplied as name-value
% pairs. Default values are in square brackets.
%
% experiment ['BinaryMixtures]: One of 'BinaryMixtures' or
% 'ComplexMixtures'; which dataset to use.
%
% cells ['PNs']: One of 'PNs' or 'KCs', which cells to use. 
%
% startingCell [1]: Integer, the ID of the cell to start with.
%
% startTime [-1]: The starting time in seconds of the raster window,
% relative to odor onset.
%
% endTime [3]: The ending time in seconds of the raster window,
% relative to odor onset.
%
% binSize [0.050]: The bin size in seconds to use in the binary
% mixtures response reconstructions (Figure 2C).
%
% showRecs [false]: Whether to show the reconstructions. This can be a
% little time-consuming, so is off by default.
%
% EXAMPLES:
%
% Plot the PNs in the binary mixtures experiments:
% >> SpikeTimesBrowser();
%
% Turn on the reconstructions:
% >> SpikeTimesBrowser('showRecs', true);
%
% Plot the results in Figure 2A-D
% >> SpikeTimesBrowser('showRecs', true, 'startingCell', 128);
% 
% Plot the PNs in the complex mixtures experiments:
% >> SpikeTimesBrowser('experiment', 'ComplexMixtures');
%
% Plot the KCs in the complex mixtures experiments:
% >> SpikeTimesBrowser('experiment', 'ComplexMixtures', 'cells', 'KCs');
%
% Look at the response over a longer time window:
% >> SpikeTimesBrowser('experiment', 'ComplexMixtures', 'cells', 'KCs', 'endTime', 5);

p = inputParser;
p.addOptional('experiment' , 'BinaryMixtures');
p.addOptional('cells',       'PNs');
p.addOptional('startingCell', 1);
p.addOptional('startTime', -1);
p.addOptional('endTime',    3);
p.addOptional('binSize', 0.050);
p.addOptional('showRecs', false);
p.parse(varargin{:});

experiment   = p.Results.experiment;
cells        = p.Results.cells;
startingCell = p.Results.startingCell;
startTime    = p.Results.startTime;
endTime      = p.Results.endTime;
binSize      = p.Results.binSize;
showRecs     = p.Results.showRecs;

odorOnsetTime = 2;

%% Validate the inputs
if (startTime >= endTime)
  error('Start time must be less than end time.');
end

if (binSize <= 0)
  error('Bin size must be greater than 0.');
end

switch lower(experiment)
 case 'binarymixtures'
  experiment = 'BinaryMixtures';
 case 'complexmixtures'
  experiment = 'ComplexMixtures';
 otherwise
  error('Experiment must be either "BinaryMixtures" or "ComplexMixtures.');
end

switch lower(cells)
 case 'pns'
  cells = 'PN';
 case 'kcs'
  if isequal(experiment, 'BinaryMixtures')
    error('No KCs recorded in the binary mixtures experiments.');
  end
  cells = 'KC';
 otherwise
  error('Cells must be either "PNs" or "KCs".');
end

if isequal(experiment, 'BinaryMixtures')
  spt       = LoadTocSpikeTimes('rawpn_binary_mixtures');
  numCells  = 168;
  numTrials = 10;
  numOdors  = 27;
  odorDur   = 0.3;
  figureTitle = 'PNs (binary mixtures)';
else
  numOdors  = 44;
  numTrials = 7;
  odorDur   = 0.5;
  if isequal(cells, 'PN')
    spt = LoadTocSpikeTimes('rawpn');
    figureTitle = 'PNs (complex mixtures)';
    numCells = 174;
  else
    spt = LoadTocSpikeTimes('rawkc');
    figureTitle = 'KCs (complex mixtures)';
    numCells = 209;
  end
end

if (startingCell<1 | startingCell > numCells)
  error('Staring cell for this dataset must be in the range 1 - %d.\n', numCells);
end

% Now make the plots
currentCell = startingCell
if (isequal(experiment, 'BinaryMixtures'))
  spt      = ConvertSpikeTimesFromSparseToFull(spt);
  PlotBinaryMixturesOdorResponseRastersForCell(currentCell, spt, startTime+odorOnsetTime, endTime+odorOnsetTime, binSize, showRecs, [])
  figureId = gcf;
  userData = struct;
  userData.currentCell = currentCell;
  userData.plotFunction = @(whichCell)  PlotBinaryMixturesOdorResponseRastersForCell(whichCell, spt, startTime+odorOnsetTime, endTime+odorOnsetTime, binSize, showRecs, figureId);
  set(figureId, 'UserData', userData,'KeyPressFcn', @BinaryMixturesKeyPressFcn);
elseif (isequal(experiment, 'ComplexMixtures'))

  if (isequal(cells, 'PN'))
    spt      = ConvertSpikeTimesFromSparseToFull(spt);
    sptCbot  = reshape(spt,[], numTrials, numOdors, numCells);    
  elseif (isequal(cells, 'KC'))
    spt      = ConvertSpikeTimesFromSparseToFull(spt);
    sptCbot  = reshape(spt,[], numTrials, numOdors, numCells);        
  else
    error('Unknown cell type "%s".\n', cells);
  end
  
  figureName = sprintf('Complex Mixtures: %s %d', cells, currentCell);
  figureId = sfigure(FindFigureCreate(figureName));
  PlotOdorResponseRastersForCell(spt, currentCell, [startTime endTime]+2, 'F', 'figureId', figureId);
  set(figureId, 'NumberTitle','off','Name', figureName);

  userData = struct;
  userData.currentCell  = currentCell;
  userData.numCells     = numCells;
  userData.plotFunction = @(whichCell) ComplexMixturesPlotFunction(cells, spt, whichCell, startTime, endTime, figureId);
  set(figureId, 'UserData', userData, 'KeyPressFcn', @ComplexMixturesKeyPressFcn);
else
  error('Unknown experiment "%s".', experiment);
end

function BinaryMixturesPlotFunction(whichCell, spt, startTime, endTime, binSize, showRecs, figureId)
PlotBinaryMixturesOdorResponseRastersForCell(currentCell, spt, startTime+odorOnsetTime, endTime+odorOnsetTime, binSize, showRecs, figureId);
set(figureId, 'name', sprintf('Binary Mixtures: PN %d', whichCell), 'NumberTitle', 'off');

function ComplexMixturesPlotFunction(cells, spt, whichCell, startTime, endTime, figureId)
PlotOdorResponseRastersForCell(spt, whichCell, [startTime endTime]+2, 'F', 'figureId', figureId);
set(figureId, 'name', sprintf('Complex Mixtures: %s %d', cells, whichCell));

function ComplexMixturesKeyPressFcn(obj, evt)
userData = get(obj,'UserData');
currentCell = userData.currentCell;
numCells    = userData.numCells;
switch evt.Key
 case 'leftarrow'
  currentCell = currentCell - 1;
  if (currentCell < 1)
    currentCell = numCells;
  end
 case 'rightarrow'
  currentCell = currentCell + 1;
  if (currentCell > numCells)
    currentCell = 1;
  end  
 case 'return'
end
userData.currentCell = currentCell;
userData.plotFunction(currentCell);
set(obj, 'UserData', userData);

function BinaryMixturesKeyPressFcn(obj, evt)
userData = get(obj,'UserData');
currentCell = userData.currentCell;
switch evt.Key
 case 'leftarrow'
  currentCell = currentCell - 1;
  if (currentCell < 1)
    currentCell = 168;
  end
 case 'rightarrow'
  currentCell = currentCell + 1;
  if (currentCell > 168)
    currentCell = 1;
  end  
 case 'return'
end
userData.currentCell = currentCell;
userData.plotFunction(currentCell);
set(obj, 'UserData', userData);
