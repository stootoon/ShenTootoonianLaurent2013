function m = TrajectoryComparisonMetrics(Xcbo,XrecCbo,whichMetric,summaryFunction)
% m = TrajectoryComparisonMetrics(Xcbo,XrecCbo,whichMetric,summaryFunction)
%
% This function is called by MAKEFIGURES and computes metrics for
% comparing trajectories to their reconstructions.
%
% whichMetric can be one of
%
% 'SSE': The sum of the squared error between the population response
% and its reconstruction,
%
% 'R2': 1 minus the ratio of SSE/SST, where SST is the sum of squares
% of the population response, at each time bin,
%
% 'CorrDist': The correlation distance between the population response
% and its reconstruction at each time bin,
%
% 'EucDist': The euclidean distance between the population response
% and its reconstruction at each time bin.
%
% Each of these functions produces a NUMBINS x NUMODORS matrix of
% metrics, computed for each of the bins for each
% odor. SUMMARYFUNCTION is applied to this matrix to produce a final
% summary metric. A typical value for it is 
%
% @(M) mean(min(M,[],1)), 
% 
% which computes the minimum for each odor and then takes the mean.

if (length(size(Xcbo))~=3)
  error('Xcbo should have three dimensions.');
end

if (size(Xcbo)~=size(XrecCbo))
  error('Xcbo and XrecCbo should have the same size.');
end

[numCells,numBins,numOdors] = size(Xcbo);

switch (lower(whichMetric))
 case 'sse'
  M = reshape(sum((Xcbo - XrecCbo).^2,1),numBins, numOdors);
 case 'r2'
  sst = reshape(sum(bsxfun(@minus, Xcbo, mean(Xcbo)).^2,1), numBins, numOdors);
  sse = reshape(sum(bsxfun(@minus, Xcbo, XrecCbo).^2,1), numBins, numOdors);
  M = 1 - sse./sst;
 case 'corrdist'
  M = corrdist(reshape(Xcbo,numCells,[]), reshape(XrecCbo,numCells,[]));
  M = reshape(M,numBins,numOdors);
 case 'eucdist'  
  M = reshape(sqrt(sum((Xcbo - XrecCbo).^2,1)), numBins, numOdors);
 otherwise
  error('Unknown metric "%s".\n', whichMetric);
end

m = summaryFunction(M);

if (length(m)~=1)
  error('Summary function should return a scalar.');
end