function Isample = ComputeTrialSamplingFromParametersStructure(Parameters)
% Isample = ComputeTrialSamplingFromParametersStructure(Parameters)
%
% Computes a cell array of vectors contains the trials to sample based
% on the contents of the PARAMETERS structure. The computation is as follows:
%
% If the IBS field Of the PARAMETERS structure is empty, it's assumed
% that no bootstrapping is desired and single trials are
% required. Hence ISAMPLE is a cell array whose elements are
% individual trials. 
% 
% If the IBS field of the PARAMETERS structure is NOT empty,
% bootstrapping is assumed and the ISAMPLES array contains the trials
% vectors contained in the columsn of the IBS field, and additionally
% the set of trial samples required for the jackknife i.e. all trials
% but the i'th, for i = 1:numTrials.

if (isempty(Parameters.Ibs))
  Isample = arrayfun(@(i) i, 1:Parameters.numTrials, 'UniformOutput', false);
else
  Ibs = mat2cell(Parameters.Ibs, size(Parameters.Ibs,1), ones(1, size(Parameters.Ibs,2)));
  Ijk = arrayfun(@(i) [1:i-1 i+1:Parameters.numTrials]', 1:Parameters.numTrials, 'UniformOutput', false);
  Isample = [Ibs Ijk];
end


