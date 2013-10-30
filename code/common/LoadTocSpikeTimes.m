function X = LoadTocSpikeTimes(which)
% X = LoadTocSpikeTimes(which)
%
% Returns the indicated matrix of Toc spike times. 
%
% Possible (case-insensitive) values of 'which' are:
%
% rawpn: The raw PNs. 7 trials, 44 odors, 174 cells.
% rawkc: The raw KCs. 7 trials, 44 odors, 209 cells.
% rawpn_binary_mixtures: The raw PNs in the binary mixtures experiments. 10 trials, 27 odors, 168 cells.

switch(lower(which))
 case 'rawpn'
  X = LoadVarFromMatFileByName(fullfile(GetRootDir('spiketimes'),'tocRawPnSpikeTimes.mat'),'tocSpikeTimes');  
  X = RemoveDuplicateCells(X,[7 44 175]);
 case 'rawpn_binary_mixtures'
  X = LoadVarFromMatFileByName(fullfile(GetRootDir('spiketimes'),'tocRawPnBinaryMixturesSpikeTimes'),'tocSpikeTimes');
 case 'rawkc'
  X = LoadVarFromMatFileByName(fullfile(GetRootDir('spiketimes'),'tocRawKcSpikeTimes.mat'),'tocSpikeTimes');  
 otherwise
  error('Did not understand input "%s".\n', which);  
end
  

