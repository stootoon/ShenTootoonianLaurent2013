function M = ComputeMetricsForTransientProximity_noduplicates(t0, binSize, responseWindow, whichComparisons, numBs)
% M = ComputeMetricsForTransientProximity_noduplicates(t0, binSize, responseWindow, whichComparisons, numBs)
%
% Computes the distances for the comparisons specified by
% WHICHCOMPARISONS for the time bins in RESPONSEWINDOW. These can be
% used to demonstrate transient proximity between odor trajectories.

t1 = CheckStartTimeAndSetEndTimeUsingResponseWindow(t0, responseWindow, binSize);
Parameters = ComputeInputParametersForTrajectoryComputations_noduplicates('PnComplexMixture','startTime', t0, 'endTime', t1, 'binSize', 0.1, 'sampleType', numBs);

t = Parameters.startTime + (0:Parameters.numBins-1)*Parameters.binSize;

whichBins = find(t>=responseWindow(1) & t<=responseWindow(2));
whichOdors = [4:11 13:44]; % Skip *high, Paraffin oil.

odorNames = GetOdorsList;
odorNames = odorNames(whichOdors);

numBins = numel(whichBins);
numCmps = numel(whichComparisons);

D = zeros(numBs, numBins, numCmps);

Isample = ComputeTrialSamplingFromParametersStructure(Parameters);
numSamples = numel(Isample);

M = struct;
M.cmps = struct;
for i = 1:numSamples
  X = squeeze(mean(Parameters.spikeCounts(:,whichBins,:,Isample{i}),4));
  for j = 1:numCmps
    [od1, od2] = PairwiseInd2Sub(numel(whichOdors), whichComparisons(j));
    M.cmps(j).odor1 = odorNames{od1};
    M.cmps(j).odor2 = odorNames{od2};
    for k = 1:numBins
      D(i,k,j) = pdist([X(:,k,od1)  X(:,k,od2)]','correlation');
    end
  end
end

M.t0 = t0;
M.firstBinTime = t(whichBins(1));
M.binSize = binSize;
M.responseWindow = responseWindow;
M.whichComparisons = whichComparisons;
M.D = D;

      
