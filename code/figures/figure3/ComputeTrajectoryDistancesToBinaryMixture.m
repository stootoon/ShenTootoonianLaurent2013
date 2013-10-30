function [Doverall, DperTimeBin] = ComputeTrajectoryDistancesToBinaryMixture(Parameters, whichMixture, whichBinsForOverall, whichBinsForPerBin, varargin)
% [Doverall, DperTimeBin] = ComputeTrajectoryDistancesToBinaryMixture(Parameters, whichMixture, whichBinsForOverall, whichBinsForPerBin, varargin) 
%
% Computes the distances from each of the mixture trajectories, whose
% data is in the structure PARAMETERS to that WHICHMIXTURE.
%
% Two types of distances are computed. DOVERALL is the distance
% between each trajectory as a whole (i.e. temporally concatenated),
% during the bins indicated in WHICBINSFOROVERALL. DPERTIMEBIN is the distances
% between the PN vectors during each of the bins in WHICHBINSFORPERBIN.
%
% PARAMETERS is a structure returned by
% COMPUTEINPUTPARAMETERSFORTRAJECTORYCOMPUTATIONS. WHICHMIXTURE is two
% element vector indicated the desired target concentration to compare
% to. WHICHBINS is a vector of desired time bins for the trajectory
% computations.
%
% DOVERALL is a NUMMIXTURES x NUMBOOTSTRAPS matrix of distances where
% DOVERALL(I,J) is the distance to the target trajectory for binary
% mixture I and bootstrap run J.
%
% DPERTIMEBIN is a NUMMIXTURES x NUMBOOTSTRAPS x NUMTIMEBINS matrix of
% distances where DOVERALL(I,J,K) is the distance to the target
% trajectory for binary mixture I and bootstrap run J, in time bin K.
%
% Distances are euclidean distances by default, but a different
% distance function, as understood by PDIST, can be supplied by the
% optional argument.

if (isempty(varargin))
  distFun = 'euclidean';
else
  distFun = varargin{1};
end

concs = [0   30  60  80  100 120 140 140 140 140 140 140 140; 
         140 140 140 140 140 140 140 120 100 80  60  30  0];

targetConcInd = GetIndexForBinaryMixtureConcentrationPair(whichMixture(1), whichMixture(2));
concInds = arrayfun(@(u,v) GetIndexForBinaryMixtureConcentrationPair(u,v), concs(1,:)', concs(2,:)');

s = Parameters.spikeCounts; % c b o t
numTrials = Parameters.numTrials; %

Isample = ComputeTrialSamplingFromParametersStructure(Parameters);
numSamples = numel(Isample);

% Now compute!
Doverall = bsxfun(@(Iodor,j) arrayfun(@(iodor) ...
                                  pdist([Columnize(mean(s(:,whichBinsForOverall,iodor,Isample{j}),4))...
                                        Columnize(mean(s(:,whichBinsForOverall,targetConcInd,  Isample{j}),4))]', distFun),...
                                      Iodor),...
                  concInds(:),1:numSamples);

DperTimeBin = cell2mat(arrayfun(@(ibin)...
                                bsxfun(@(Iodor,isamp)...
                                       arrayfun(@(iodor)...
                                                pdist(...
                                                    [squeeze(mean(s(:,ibin,iodor,Isample{isamp}),4)) ...
                                                     squeeze(mean(s(:,ibin,targetConcInd,Isample{isamp}),4))]',... 
                                                    distFun),...
                                                Iodor),...
                                       concInds(:), 1:numSamples),...
                                whichBinsForPerBin,'UniformOutput',false));

DperTimeBin = reshape(DperTimeBin, numel(concInds), numSamples, numel(whichBinsForPerBin));
DperTimeBin = permute(DperTimeBin, [1 3 2]); % bootstrap index should be last.
                                       