function [Y, binStarts] = CountSpikesInBinsAndAverageAcrossTrials(X, whichTrials, whichOdors, whichCells, varargin)
% [Y, binStarts] = CountSpikesInBinsAndAverageAcrossTrials(X, whichTrials, whichOdors, whichCells, varargin)
%
% Reads the TOC matrix X of spike times and extracts the subset
% specified in the arrays WHICHTRIALS, WHICHODORS, and WHICHCELLS,
% counts spikes in the binned odor responses, averages across the
% specified trials, and returns a numCells x numBins x numOdors x
% numAverages matrix Y with the result.
%
% Y(i,j,k,l) contains the average number of spikes from cell i in bin
% j to odor k for average l.
%
% If whichTrials is a cell array, each element is interpreted as a set
% of trials to average over.
%
% Optional arguments can be provided in name-value pairs:
%
% startTime: First bin time in seconds. Default is 0.
% endTime: Time of end of last bin in seconds. Default is 2.
% binSize: Bin size in seconds. Default is 0.050.
% numAllTrials: Total number of trials. Default is 7.
% numAllOdors: Total number of odors. Default is 44.
%
% Example: 
%
% X = LoadTocSpikeTimes('rawpn');
%
% % Grab the first three odors
% whichOdorInds = 1:3;
% 
% % Average over the first and second trials.
% whichTrials = [1 2];
%
% % Use all cell.s
% whichCells = 1:174;
%
% Y = CountSpikesInBinsAndAverageAcrossTrials(X, whichTrials, whichOdorInds, whichCells);
% Y = reshape(Y, 174, []);
% L = lle(Y,K,dmax);
% L = reshape(L,[],numBins, numOdors);  
% sfigure(1); clf;
% PlotDimensionReductionResults(L,whichOdors,doSmooth);
Options = LoadStructFromNameValuePairs(varargin,...
                                       {'startTime','endTime', 'binSize','binOverlap','numAllTrials', 'numAllOdors'},...
                                       {0, 2, 0.050, 0, 7, 44});

UnpackStructureFieldsAsVariables(Options); % The fields in Options are now variables.

% Convert X to a full matrix 
X = full(X); X(X == 0) = -Inf;

numAllCells = size(X,2)/numAllTrials/numAllOdors;

% Flatten whichTrials for validation.
if (iscell(whichTrials)) 
  allTrials = Columnize(cell2mat(whichTrials));
else
  allTrials = Columnize(whichTrials);
end

%% Validate Inputs
IP = inputParser;

IP.addRequired('X',            @(x) validateattributes(x, {'numeric'}, {'2d','nonnan','nonempty'}));
IP.addRequired('whichTrials',  @(x) validateTrials(x, numAllTrials));
IP.addRequired('whichOdors',   @(x) validateattributes(x, {'numeric'}, {'<=',numAllOdors, '>=', 0}));
IP.addRequired('whichCells',   @(x) validateattributes(x, {'numeric'}, {'<=',numAllCells, '>=', 0}));

IP.addParamValue('startTime',    0,     @(x) validateattributes(x, {'numeric'}, {'scalar','finite','nonnan'}));
IP.addParamValue('endTime',      2,     @(x) validateattributes(x, {'numeric'}, {'scalar','finite','nonnan','>=', startTime}));
IP.addParamValue('binSize',      0.050, @(x) validateattributes(x, {'numeric'}, {'scalar','finite','nonnan','positive'}));
IP.addParamValue('numAllTrials', 7,     @(x) validateattributes(x, {'numeric'}, {'scalar','finite','nonnan','positive'}));
IP.addParamValue('numAllOdors',  44,    @(x) validateattributes(x, {'numeric'}, {'scalar','finite','nonnan','positive'}));

IP.parse(X, whichTrials, whichOdors, whichCells, varargin{:});

%% Inputs are validated, start the computation.
numOdors  = numel(whichOdors);
numCells  = numel(whichCells);

binStarts = startTime:binSize:endTime+binSize; % +binSize so that below when we look for the lastBin, we always find it.

% We compute binEnds the way we do, as a staggered version of binStarts 
% instead of binStarts + binSize, because it reduces roundoff errors which lead
% errorneous spike counts. For example, if we have a spike at exactly 2.3,
% if there's roundoff error, the bound of the semiopen intervals will be slightly
% greater than 2.3, and we'll count a spike in [2.2, 2.3), when we shouldn't.
% 
% But if we use one set of values to mark the bin boundaries, then we
% can at least avoid problems of counting spikes twice because we avoid overlap
% between bins due to roundoff errors.
%
% Due to this change, we've stopped allowing overlap between bins, otherwise it 
% gets a bit complicated choosing where bins start and end.

binEnds   = binStarts(2:end);

lastBin   = find(binEnds >= endTime-10*eps);
binStarts = binStarts(1:lastBin);
binEnds   = binEnds(1:lastBin);
numBins   = numel(binStarts);

if (~iscell(whichTrials))
  whichTrials = {whichTrials};
end

numAverages = numel(whichTrials);

Y = zeros(numCells, numOdors, numBins, numAverages);
Y = reshape(Y, numCells, []);
ip = 1; % insertion point

for i = 1:numAverages
  trials = whichTrials{i};
  Xsubset = GetSubsetOfTocMatrix(X, [numAllTrials numAllOdors numAllCells], ...
                                 {trials, whichOdors, whichCells});
  for j = 1:numBins % for each time bin
    t0 = binStarts(j);
    t1 = binEnds(j);
    C = CountSpikesInSemiClosedTimeWindow(Xsubset,t0,t1); % Count number of spikes in each trial
    C = reshape(C, numel(trials), numOdors*numCells);
    C = mean(C,1); % Average across trials
    C = reshape(C, numOdors, numCells)';
    Y(:,ip+(0:numOdors-1)) = C; % Stick in matrix
    ip = ip + numOdors;
  end

end
Y = reshape(Y,numCells, numOdors, numBins, numAverages);
Y = permute(Y,[1 3 2 4]);

function valid = validateTrials(whichTrials, maxTrial)
valid = 0;

if (iscell(whichTrials))
  if (any(cellfun(@iscell, whichTrials))); % expecting numeric array, not more cells.
    error('Trials list must contain only numeric arrays.');
  elseif (any(cellfun(@isempty, whichTrials)))
    error('Some trials lists are empty.');
  else
    whichTrials = cell2mat(whichTrials);
  end
end

assert(~isempty(whichTrials),'Trials list is empty.');
assert(~any(isnan(whichTrials) | isinf(whichTrials)), 'Some trial values are NaN or Inf.');
assert(~any(whichTrials < 1 | whichTrials>maxTrial),  'Some trial values are out of range.');
valid = 1;
