function [mixtureInds, mixtureVals] = GetBinaryMixtureIndsForResponseLinearity()
% function [mixtureInds, mixtureVals] = GetBinaryMixtureIndsForResponseLinearity()
%
% Returns the mixture indices and mixture values used in the response
% linearity fits i.e. the binary mixture morphs and the concentration
% series.
[mixtureInds, mixtureVals] = GetBinaryMixtureMorphInds;
mixtureVals = mixtureVals(2:end-1,:); % Use only the ones that have some of each.
mixtureInds = mixtureInds(2:end-1);
numMorphOdors = length(mixtureInds);
concentrationSeries = [30 60 80 100 140]'*[1 1];
for i = 1:size(concentrationSeries,1)
  mixtureVals(numMorphOdors+i,:) = concentrationSeries(i,:);
  mixtureInds(numMorphOdors+i)  = GetIndexForBinaryMixtureConcentrationPair(concentrationSeries(i,1), concentrationSeries(i,2));
end
