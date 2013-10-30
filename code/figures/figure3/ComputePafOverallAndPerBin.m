function [PAFov, PAFpb] = ComputePafOverallAndPerBin(Parameters, whichBins, mixtureInds)
% function [PAFov, PAFpb] = ComputePafOverallAndPerBin(Parameters, whichBins, mixtureInds)
%
% Computes the PAFS relative to the first mixture of the thirteen
% whose indices should be in MIXTUREINDS. The projection angle
% fractions are computed in two ways: globally and locally. For the
% global computation, the spike counts for all the time bins specified
% in WHICHBINS are concatenated to yield one vector per mixture, for
% which the PAFS are then computed and returned in the NUMMIXTURES-2 x
% NUMBOOTSTRAPS matrix PAFOV. The PAFs are also computed locally, one
% bin at a time for each of the bins in WHICHBINS, and the results
% returned in the NUMMIXTURES-2 x NUMBINS x NUMBOOTSTRAPS matrix PAFPB.
%
% PAFS are computed relative to the first mixture in MIXTUREINDS.
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
PAFpb = computePAF(U0,U2,nonPureInd, @(i) reshape(Y(:,:,i,:),size(Y,1),[]));
PAFpb = reshape(PAFpb, numNonPure, numBins, numSamples);

%% PAFS GLOBALLY
Xov = reshape(Y, [], numMix, numSamples); % Concatenate the response window bins.
U0 = reshape(Xov(:,firstCmpInd,:), size(Xov,1),[]);
U2 = reshape(Xov(:,secondCmpInd,:), size(Xov,1),[]);
PAFov = computePAF(U0,U2,nonPureInd, @(i) reshape(Xov(:,i,:), size(Xov,1),[]));
PAFov = reshape(PAFov, numNonPure, numSamples);

function PAF = computePAF(U0,U2,mixtureInds,f_getMixture)
% Helper functions
PROJ       = @(u,v,w) [v w]*([v w]\u); % Returns the projection of U in the span of v and w.
PROJ_ANGLE = @(U,v)   arrayfun(@(i) acos(dot(U(:,i),v(:))/norm(U(:,i))/norm(v(:))),1:size(U,2)); % Computes the angle between the each column of U into the span of v.
NORM       = @(U)     sqrt(sum(U.^2,1)); % The Euclidean norm of each column.

% Compute the intercomponent angles
N0 = NORM(U0);
N2 = NORM(U2);
A0 = acos(sum(U0.*U2,1)./N0./N2);

% Compute the projection angles
A1 = zeros(numel(mixtureInds),size(A0,2));
for i = 1:numel(mixtureInds)
  X1 = f_getMixture(mixtureInds(i));
  % Compute the projections
  U1 = cell2mat(arrayfun(@(j) PROJ(X1(:,j), U0(:,j), U2(:,j)),1:size(X1,2),'UniformOutput',false));
  N1 = NORM(U1);
  % Compute the projection angle
  A1(i,:) = acos(sum(U1.*U0,1)./NORM(U1)./NORM(U0));
end
PAF = bsxfun(@rdivide, A1, A0);
