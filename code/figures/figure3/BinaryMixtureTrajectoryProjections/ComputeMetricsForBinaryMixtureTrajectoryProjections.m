function M = ComputeMetricsForBinaryMixtureTrajectoryProjections(t0, binSize, responseWindow, numBs)
% function M = ComputeMetricsForBinaryMixtureTrajectoryProjections(t0, binSize, responseWindow, numBs)
%
% For the given arguments, creates the PARAMETERS structure and
% computes the following mixture morph metrics to address whether the
% extent to which the trajectory for the 1:1 mixture of octanol and
% citral lies in between the component trajectories:
%
% 1. Global distance to citral (returned in M.DOVCIT);
% 2. Perbin distance to citral (M.DPBCIT)
% 3. Global PAF to citral (M.PAFOVCIT)
% 4. Perbin PAF to citral (M.PAFPBCIT)
% 5. Global distance to octanol (returned in M.DOVOCT);
% 6. Perbin distance to octanol (M.DPBOCT)
% 7. Global PAF to octanol (M.PAFOVOCT)
% 8. Perbin PAF to octanol (M.PAFPBOCT)
% 9. Global PMF (M.PMFOV)
% 10. Perbin PMF (M.PMFPB)
%
% M additional contains fields containing the provided input
% arguments, as they were used to compute the Parameters structure for
% which the metrics were computed, and also the IBS sub field of
% parameters indicating the trials used in forming the bootstrap
% trajectories.

if (t0>responseWindow(1))
  error('Start time (T0) is greater than start of response window.');
end
t1 = responseWindow(2)+2*binSize;

Parameters = ComputeInputParametersForTrajectoryComputations_noduplicates('BinaryMixture','startTime',t0, 'endTime',t1, 'binSize', binSize,'sampleType', numBs);

P = Parameters;
t = P.startTime+(0:P.numBins-1)*P.binSize;

whichBins = find(t>=responseWindow(1) & t<responseWindow(2));

mixtureVals = [0 140; 30 140; 60 140; 80 140; 100 140; 120 140; 140 140;
               140 120; 140 100; 140 80; 140 60; 140 30; 140 0];
mixtureInds = arrayfun(@(i) GetIndexForBinaryMixtureConcentrationPair(mixtureVals(i,1),mixtureVals(i,2)),1:13);

M = struct;
[M.PAFovCit, M.PAFpbCit] = ComputePafOverallAndPerBin(P, whichBins, mixtureInds);
[M.PAFovOct, M.PAFpbOct] = ComputePafOverallAndPerBin(P, whichBins, fliplr(mixtureInds));
[M.DovCit,   M.DpbCit]   = ComputeTrajectoryDistancesToBinaryMixture(P, [0 140], whichBins, whichBins);
[M.DovOct,   M.DpbOct]   = ComputeTrajectoryDistancesToBinaryMixture(P, [140 0], whichBins, whichBins);
[M.PMFov, M.PMFpb]       = ComputePmfOverallAndPerBin(P, whichBins, mixtureInds);
M.t0 = t0;
M.binSize = binSize;
M.responseWindow = responseWindow;
M.Ibs = Parameters.Ibs;
