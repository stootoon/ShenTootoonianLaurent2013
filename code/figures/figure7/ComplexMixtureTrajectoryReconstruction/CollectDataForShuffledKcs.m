function [Vsh,VshCbo, Ush, UshCbo, Bsh, b0sh] = CollectDataForShuffledKcs(Results)
% [Vsh VshCbo, Ush, UshCbo, Bsh, b0sh] = CollectDataForShuffledKcs(Results)
%
% Given the cell array of shuffles RESULTS, collects the results and
% returns:
%
% Vsh: A NUMSAMPLES x NUMKCS x NUMSHUFFLES matrix of shuffled KC responses,
% VshCbo: Vsh in NUMKCS x NUMBINS x NUMODORS x NUMSHUFFLES format.
% Ush: A NUMSAMPLES x NUMPNS x NUMSHUFFLES matrix of PN reconstructions.
% UshCbo: Ush in NUMPNS x NUMBINS x NUMODORS x NUMSHUFFLES format.
% Bsh: A NUMKCS x NUMPNS x NUMSHUFFLES matrix of reconstruction weights,
% b0sh: A NUMSHUFFLES X NUMPNS matrix of reconstruction offsets.

numShuffles = length(Results);

Vsh = repmat(Results{1}{1},[1,1,numShuffles]);
Bsh = repmat(Results{1}{2},[1,1,numShuffles]);
b0sh= repmat(Results{1}{3},numShuffles,1);
[numSamples, numKcs] = size(Results{1}{1});
[numKcs, numPns] = size(Results{1}{2});

numOdors = length(GetOdorsList);
numBinsPerOdor = numSamples/numOdors;

VshCbo = zeros(numKcs, numBinsPerOdor, numOdors, numShuffles);

Ush = zeros(numSamples, numPns, numShuffles);
UshCbo = zeros(numPns, numBinsPerOdor, numOdors, numShuffles);

for i = 1:numShuffles
  Vsh(:,:,i) = Results{i}{1};
  VshCbo(:,:,:,i) = permute(reshape(Vsh(:,:,i),[numBinsPerOdor,numOdors,numKcs]),[3 1 2]);
  Bsh(:,:,i) = Results{i}{2};
  b0sh(i,:) = Results{i}{3};
  Ush(:,:,i) = bsxfun(@plus,Vsh(:,:,i)*Bsh(:,:,i),b0sh(i,:));
  UshCbo(:,:,:,i) = permute(reshape(Ush(:,:,i),[numBinsPerOdor, numOdors, numPns]),[3 1 2]);
end


