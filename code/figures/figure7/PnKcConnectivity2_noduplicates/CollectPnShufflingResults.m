function [Ush,Bsh,b0sh] = CollectPnShufflingResults(Results)
% [Ush,Bsh,b0sh] = CollectPnShufflingResults(Results)
%
% Given the cell array of PnShuffling results RESULTS (as produced by
% PARCOMPUTELASSOWEIGHTSFORSHUFFLED, collects them into the 
%
% NUMSAMPLES x NUMPNS x NUMSHUFFLES matrix Ush, 
% NUMPNS x NUMKCS x NUMSHUFFLES weights matrix Bsh,
% NUMSHUFFLES x NUMKCS offsets matrix b0sh.

Results = Results(:);
numShuffles = numel(Results);

R = Results{1};
[numSamples, numPns] = size(R{1});
[numPns, numKcs] = size(R{2});

Ush = zeros(numSamples, numPns, numShuffles);
Bsh = zeros(numPns, numKcs, numShuffles);
b0sh= zeros(numShuffles, numKcs);

for i = 1:numShuffles
  Ush(:,:,i) = Results{i}{1};
  Bsh(:,:,i) = Results{i}{2};
  b0sh(i,:)  = Results{i}{3};
end