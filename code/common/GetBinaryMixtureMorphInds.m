function [mixtureInds, mixtureVals] = GetBinaryMixtureMorphInds()
% [mixtureInds, mixtureVals] = GetBinaryMixtureMorphInds()
%
% Returns the indices (in the full list of mixtures used in the binary
% mixtures experiments) of the mixture pairs used in the mixture morph
% experiments.

% The mixture values used in the morph experiments
mixtureVals = [0 140; 30 140; 60 140; 80 140; 100 140; 120 140; 140 140;
               140 120; 140 100; 140 80; 140 60; 140 30; 140 0;];

mixtureInds = arrayfun(@(i) GetIndexForBinaryMixtureConcentrationPair(mixtureVals(i,1), mixtureVals(i,2)),1:13);
