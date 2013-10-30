function [Vsh,Bsh,b0sh] = CollectKcShufflingResults(Results)
% [Vsh,Bsh,b0sh] = CollectKcShufflingResults(Results)
%
% Given the cell array of KcShuffling results RESULTS (as produced by
% PARCOMPUTELASSOWEIGHTSFORSHUFFLEDKCS, collects them into the 
%
% NUMSAMPLES x NUMKCS x NUMSHUFFLES matrix Vsh, 
% NUMPNS x NUMKCS x NUMSHUFFLES weights matrix Bsh,
% NUMSHUFFLES x NUMKCS offsets matrix b0sh.

Results = Results(:);
numShuffles = numel(Results);

R = Results{1};
[numSamples, numPns] = size(R{1});
[numPns, numKcs] = size(R{2});

Vsh = zeros(numSamples, numKcs, numShuffles);
Bsh = zeros(numPns, numKcs, numShuffles);
b0sh= zeros(numShuffles, numKcs);

for i = 1:numShuffles
  Vsh(:,:,i) = Results{i}{1};
  Bsh(:,:,i) = Results{i}{2};
  b0sh(i,:)  = Results{i}{3};
end