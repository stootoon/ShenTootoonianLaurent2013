function ind =  GetResponseLinearityIndexForBinaryMixtureConcentrationPair(whichOct, whichCit)
% ind =  GetResponseLinearityIndexForBinaryMixtureConcentrationPair(whichOct, whichCit)
%
% Given the concentrations of citral and octanol, returns the index of
% the corresponding mixture as used in the Binary Mixture Response
% linearity analyses. The difference between this indexing and the
% standard is e.g. that pure mixtures aren't included, and the morph
% series is first, followed by the concentration series, and that
% 140:140 appears twice. An example function that uses this index is
% FITALLMODELSFORCELLANDMIXTUREBAYES1LAPLACE.

[mixtureInds, mixtureVals] = GetBinaryMixtureIndsForResponseLinearity;

ind = find(arrayfun(@(i) isequal(mixtureVals(i,:), [whichOct whichCit]), 1:numel(mixtureInds)), 1);
