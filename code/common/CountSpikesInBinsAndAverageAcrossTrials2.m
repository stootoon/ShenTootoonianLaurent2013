function [Y, binStarts, Imove] = CountSpikesInBinsAndAverageAcrossTrials2(X, whichTrials, whichOdors, whichCells, varargin)
% [Y, binStarts, Imove] = CountSpikesInBinsAndAverageAcrossTrials2(X, whichTrials, whichOdors, whichCells, varargin)
%
% Same as COUNTSPIKESINBINSANDAVERAGEACROSSTRIALS, but allows
% overlapping windows. The problem of counting spikes near bin edges
% twice is solved by moving spikes by a small amount so that they fall
% clearly within one bin or another. The overlap amount is provided as
% the optional argumetn 'overlap', in seconds.
%
% See also: CountSpikesInBinsAndAverageAcrossTrials.

Options = LoadStructFromNameValuePairs(varargin,...
                                       {'startTime','endTime', 'binSize','overlap','numAllTrials', 'numAllOdors'},...
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
IP.addParamValue('overlap',      0,     @(x) validateattributes(x, {'numeric'}, {'scalar','finite','nonnan'}));
IP.addParamValue('deltaSpike',   1e-9,  @(x) validateattributes(x, {'numeric'}, {'scalar','finite','nonnan','positive'}));
IP.addParamValue('numAllTrials', 7,     @(x) validateattributes(x, {'numeric'}, {'scalar','finite','nonnan','positive'}));
IP.addParamValue('numAllOdors',  44,    @(x) validateattributes(x, {'numeric'}, {'scalar','finite','nonnan','positive'}));

IP.parse(X, whichTrials, whichOdors, whichCells, varargin{:});

if (IP.Results.overlap>= IP.Results.binSize)
  error('Overlap must be < binsize.');
end

%% Inputs are validated, start the computation.
numOdors  = numel(whichOdors);
numCells  = numel(whichCells);

overlap = max(IP.Results.overlap,0);
if (1-overlap/binSize<1e-4)
  overlap = 0;
end
stepSize = binSize - overlap;
binStarts = startTime:stepSize:endTime-binSize;
binEnds   = min(binStarts+binSize,endTime);
numBins   = numel(binStarts);

% Now find all spikes that are near binEdges and move them
binEdges = unique([binStarts binEnds]);
deltaSpike = IP.Results.deltaSpike;

% Make those slightly ahead slightly more ahead
Imove = sparse(size(X));
for i = 1:size(X,2)
  dX = bsxfun(@minus,X(:,i),binEdges);
  
  [iahead,column] = find(dX>0 & dX<=deltaSpike);
  X(iahead,i) = X(iahead,i)+deltaSpike;
  Imove(iahead,i) = 1;

  [ibehind,column] = find(dX<0 & dX>=-deltaSpike);
  X(ibehind,i) = X(ibehind,i)-deltaSpike;
  Imove(ibehind,i) = 1;
end

% Continue with the rest of the function
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
