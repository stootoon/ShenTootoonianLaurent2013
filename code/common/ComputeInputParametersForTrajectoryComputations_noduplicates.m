function Parameters = ComputeInputParametersForTrajectoryComputations_noduplicates(dataset, varargin)
% Parameters = ComputeInputParametersForTrajectoryComputation(dataset, ...)
%
% Computes the parameters used in the computation of the bootstrap
% trajectories. DATASET specifies the dataset from which the spike
% counts are computed. It can be one of
% 'BinaryMixture','PnComplexMixture','KcComplexMixture'. The
% returned parameter structure has the following fields. Those marked
% with asterixes can be overridden by providing name-value
% pairs. Default values for these are indicated in square brackets.
%
% DATASET: The specified data set used.
% SAMPLETYPE(*): [NUMTRIALS,CEIL(NUMTRIALS/2)].  If only one argument is provide, it's
% interpreted as the number of bootstrap trajectories to generate, by
% sampling the trials with replacement. If two arguments, N and M are
% provided, then the trials are sampled without replacement, and
% indices are generated for all unique, order-insensitive subsets of M
% elements of the set 1:N. The SAMPLETYPE field is then set to
% 'Bootstrap' or 'Subset', respectively.
%
% STARTTIME(*): [2] The start of the first time bin. Odor onset is at 2.
% ENDTIME(*)  : [5] The end of the last time bin. Odor offset is at 2.5.%
% BINSIZE(*)  : [0.050] The size of the bins to use for counting binning spike counts, in seconds.
%
% IBS: A NUMTRIALS x NUMBOOTSTRAPRUNS matrix of trial indices to use
% in each bootstrap run. 
% NUMBOOTSTRAPRUNS: The number of bootstrap runs generated.
% NUMTRIALS: The number of total trials, read from the dataset.
% NUMMIXTURES: The number of mixtures, read from the dataset.
% NUMCELLS: The number of cells, read from the dataset.
% NUMBINS: The number of bins. 
% SPIKECOUNTFORMAT: The format of the binned spike counts matrix, 'cbot'.
% SPIKECOUNTS: The binned spike counts.
% CREATEDDATE: A string containing the time the structure was created.
%
%
% Example: Compute the parameters for start time 2.02, end time 5.03,
% bin size 0.045, using all unique samples of 2 of the 10 trials for
% the PN binary mixture data.
%
% >> P = ComputeInputParametersForTrajectoryComputations('BinaryMixture',
%     'sampleType',[10 2], 'startTime',2.02,'endTime',5.03,'binSize',0.045);
% 
% Generating all unique sets of indices for sampling 2 elements from the set 1:10.
% Loading spike counts for [2.020:0.045:5.030]...done.
%
% P = 
%
%        startTime: 2.0200
%          endTime: 5.0300
%          binSize: 0.0450
%       sampleType: 'Subset'
%              Ibs: [2x45 double]
% numBootstrapRuns: 45
%        numTrials: 10
%      numMixtures: 27
%         numCells: 168
% spikeCountFormat: 'cbot'
%      spikeCounts: [4-D double]
%          numBins: 67
%      createdDate: '11-Jul-2011 19:15:59'

% Load the varargs once to get what type of dataset we want.

switch(lower(dataset))
 case 'binarymixture'
  spt = ConvertSpikeTimesFromSparseToFull(LoadTocSpikeTimes('rawpn_binary_mixtures'));
  numTrials   = 10;
  numMixtures = 27;
  numCells    = 168;
 case 'pncomplexmixture'
  spt = ConvertSpikeTimesFromSparseToFull(LoadTocSpikeTimes('rawpn'));
  numTrials   = 7;
  numMixtures = 44;
  numCells    = 174;
 case 'kccomplexmixture'
  spt = ConvertSpikeTimesFromSparseToFull(LoadTocSpikeTimes('rawkc'));
  numTrials   = 7;
  numMixtures = 44;
  numCells    = 209;
 otherwise
  error('Unknown dataset: "%s".', dataset);
end

% Now read the optional parameters, setting the defaults to match 
Parameters = LoadStructFromNameValuePairs(varargin,...
                                      {'dataset', 'sampleType', 'startTime', 'endTime', 'binSize'},...
                                      {'BinaryMixture',[numTrials ceil(numTrials/2)], 2, 5, 0.050});
Parameters.dataset = lower(dataset);
switch(numel(Parameters.sampleType))
 case 1
  numBootstrapRuns = Parameters.sampleType;
  if (numBootstrapRuns == 0) % Just do the 
    Parameters.sampleType = 'Single trials';
    Ibs = [];
  else
    fprintf('Generating %d bootstrap runs by sampling trials 1:%d with replacement.\n', numBootstrapRuns, numTrials);
    ii = (1:numTrials)';
    Ibs = [ii(randi(numTrials, numTrials, numBootstrapRuns))];
    Ibs = [ii Ibs]; % Append the 'unbootstrapped' ordering at the beginning.
    Parameters.sampleType = 'Bootstrap';
  end
 case 2
  numTotal = Parameters.sampleType(1);
  numSubset = Parameters.sampleType(2);
  fprintf('Generating all unique sets of indices for sampling %d elements from the set 1:%d.\n', numSubset, numTotal);
  Ibs = GenerateUniqueResamplingWithoutReplacementIndices(numTotal, numSubset)';
  numBootstrapRuns = size(Ibs,2);
  Parameters.sampleType = 'Subset';
 otherwise
  error('Expected exactly 1 or 2 input arguments.');
end

Parameters.Ibs = Ibs;
Parameters.numBootstrapRuns = numBootstrapRuns;
Parameters.numTrials = numTrials;
Parameters.numMixtures = numMixtures;
Parameters.numCells = numCells;
Parameters.spikeCountFormat = 'cbot';
fprintf('Loading spike counts for [%1.3f:%1.3f:%1.3f]...', Parameters.startTime, Parameters.binSize, Parameters.endTime);
Parameters.spikeCounts = CountSpikesInBinsAndAverageAcrossTrials(spt, map(@Identity, 1:numTrials), 1:numMixtures, 1:numCells,'startTime',Parameters.startTime,'endTime',Parameters.endTime,'numAllTrials',numTrials,'numAllOdors',numMixtures,'binSize',Parameters.binSize);
fprintf('done.\n');
Parameters.numBins = size(Parameters.spikeCounts,2);
Parameters.createdDate = datestr(now);
