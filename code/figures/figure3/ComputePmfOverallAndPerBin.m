function [PMFov, PMFpb] = ComputePmfOverallAndPerBin(Parameters, whichBins, mixtureInds)
% function [PMFov, PMFpb] = ComputePmfOverallAndPerBin(Parameters, whichBins, mixtureInds)
%
% Computes the projection magnitude fractions for each of the mixtures
% relative to the components. The projection magnitude fractions are
% computed in two ways: globally and locally. For the global
% computation, the spike counts for all the time bins specified in
% WHICHBINS are concatenated to yield one vector per mixture, for
% which the PMFs are then computed and returned in the NUMMIXTURES-2 x
% NUMBOOTSTRAPS matrix PMFOV. The PMFs are also computed locally, one
% bin at a time for each of the bins in WHICHBINS, and the results
% returned in the NUMMIXTURES-2 x NUMBINS x NUMBOOTSTRAPS matrix
% PMFPB.
%
% The indices of the pure components are expected to be first and last
% in mixtureInds.
%
% PARAMETERS should be a structure of spike counts and other
% parameters computed by
% COMPUTEINPUTPARAMETERSFORTRAJECTORYCOMPUTATIONS. WHICHBINS is a
% vector of bins indexing the spike counts matrix in PARAMETERS
% specifying the bins to use for the local and global
% computation. MIXTUREINDS is a vector of indices for the mixtures for
% indexing the spike counts matrix in the PARAMETERS structure. It is
% expected that the first and last index are for the components, and
% those in between are for the intermediate mixtures.

firstCmpInd  = 1;
secondCmpInd = numel(mixtureInds);
nonPureInd = firstCmpInd+1:secondCmpInd-1;

Ibs        = Parameters.Ibs;
numIbs     = size(Ibs,2); % = 1 + numBs, to include the unbootstrapped ordering at the beginning.
numBins    = numel(whichBins);
numMix     = numel(mixtureInds);
numTrials  = Parameters.numTrials;
numNonPure = numel(nonPureInd);

numJackknife = numTrials;

X = Parameters.spikeCounts(:, whichBins, mixtureInds, :); % cbot format

Isamples = ComputeTrialSamplingFromParametersStructure(Parameters);
numSamples = numel(Isamples);

Y = zeros(size(X,1),size(X,2),size(X,3), numSamples);
for i = 1:numSamples
  Y(:,:,:,i) = mean(X(:,:,:,Isamples{i}),4);
end

%% PAFS PER BIN
U0 = reshape(Y(:,:,firstCmpInd,:),  size(Y,1),[]);
U2 = reshape(Y(:,:,secondCmpInd,:), size(Y,1),[]);
PMFpb = computePMF(U0,U2,nonPureInd, @(i) reshape(Y(:,:,i,:),size(Y,1),[]));
PMFpb = reshape(PMFpb, numNonPure, numBins, numSamples);

%% PMFS GLOBALLY
Xov = reshape(Y, [], numMix, numSamples); % Concatenate the response window bins.
U0 = reshape(Xov(:,firstCmpInd,:), size(Xov,1),[]);
U2 = reshape(Xov(:,secondCmpInd,:), size(Xov,1),[]);
PMFov = computePMF(U0,U2,nonPureInd, @(i) reshape(Xov(:,i,:), size(Xov,1),[]));
PMFov = reshape(PMFov, numNonPure, numSamples);

function PMF = computePMF(U0,U2,mixtureInds,f_getMixture)
% Helper functions
PROJ       = @(u,v,w) [v w]*([v w]\u); % Returns the projection of U in the span of v and w.
NORM       = @(U)     sqrt(sum(U.^2,1)); % The Euclidean norm of each column.

% Compute the projections
PMF = zeros(numel(mixtureInds),size(U0,2));
for i = 1:numel(mixtureInds)
  X1 = f_getMixture(mixtureInds(i));
  
  % Compute the projections
  X1p = cell2mat(arrayfun(@(j) PROJ(X1(:,j), U0(:,j), U2(:,j)),1:size(X1,2),'UniformOutput',false));
  N1p = NORM(X1p);
  N1  = NORM(X1);

  % Compute the projection magnitude fraction
  PMF(i,:) = N1p./N1;
end

