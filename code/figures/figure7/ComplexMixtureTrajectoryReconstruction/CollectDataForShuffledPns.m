function [Ush, UshCbo, Bsh, b0sh, UshRec, UshRecCbo] = CollectDataForShuffledPns(Results, V)
% [Ush, UshCbo, Bsh, b0sh, UshRec, UshRecCbo] = CollectDataForShuffledPns(Results, V)
%
% Given the cell array of shuffles RESULTS, and a NUMSAMPLES x NUMKCS
% matrix of KC responses, collects the results and returns:
%
% Ush: A NUMSAMPLES x NUMPNS x NUMSHUFFLES matrix of shuffled PNs
%
% UshCbo: Ush in NUMPNS x NUMBINS x NUMODORS x NUMSHUFFLES format.
%
% UshRec: A NUMSAMPLES x NUMPNS x NUMSHUFFLES matrix contain the reconstructions of the above trajectories
% UshRecCbo: UshRec in NUMPNS x NUMBINS x NUMODORS x NUMSHUFFLES format.
%
% Bsh: A NUMKCS x NUMPNS x NUMSHUFFLES matrix of reconstruction
% weights, to reconstruct the shuffled PNs with the unshuffled KCs.
%
% b0sh: A NUMSHUFFLES X NUMPNS matrix of reconstruction
% offsets.

numShuffles = length(Results);
[numSamples, numPns] = size(Results{1}{1});
[numKcs, numPns] = size(Results{1}{2});
numOdors = length(GetOdorsList);
numBinsPerOdor = numSamples/numOdors;

Ush = repmat(Results{1}{1},[1 1 numShuffles]);
UshRec = 0*Ush;
Bsh = repmat(Results{1}{2},[1 1 numShuffles]);
b0sh = repmat(Results{1}{3},[numShuffles, 1]);

for i = 1:numShuffles  
  Ush(:,:,i) = Results{i}{1};
  Bsh(:,:,i) = Results{i}{2};
  b0sh(i,:)  = Results{i}{3};
  UshRec(:,:,i) = bsxfun(@plus, V*Bsh(:,:,i), b0sh(i,:));
end

UshCbo    = permute(reshape(Ush,   numBinsPerOdor,numOdors,numPns,numShuffles),[3 1 2 4]);
UshRecCbo = permute(reshape(UshRec,numBinsPerOdor,numOdors,numPns,numShuffles),[3 1 2 4]);
