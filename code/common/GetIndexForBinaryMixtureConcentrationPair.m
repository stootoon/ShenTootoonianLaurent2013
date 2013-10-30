function ind = GetIndexForBinaryMixtureConcentrationPair(c1,c2)
% ind = GetIndexOfBinaryMixtureMixtureConcentrationPair(c1,c2)
%
% Returns an index for the specified concentration pair. This index is
% useful for accessing toc matrices of spike times. 
%
% C1 is the concentration of octanol, and C2 is the concentration of
% citrol.
%
% Example: Get the data for the mixture of 140 octanol to 80 citrol.
%
% mixtureInd = GetIndexForBinaryMixtureConcentrationPair(140,80);
% pnSpt = LoadTocSpikeTimes('rawpn_binary_mixtures');
% pnSpt = ConvertSpikeTimesFromSparseToFull(pnSpt);
% X = GetSubsetOfTocMatrix(pnSpt, [10 27 168], {[], mixtureInd, []});

pairedConcs = GetBinaryMixturePairedConcentrations;
ind = find(pairedConcs(:,1) == c1 & pairedConcs(:,2) == c2);
if (isempty(ind))
  error('Could not find concentration pair (%d,%d).', c1, c2);
end

