function Fit = FitAllModelsForCellAndMixtureBayes1Laplace(whichCell, whichMixture, varargin)
% Fit = FitAllModelsForCellAndMixture(whichCell, whichMixture, varargin)
%
% Loads the data for the specified cell and mixture, performs the fits
% for all the models, and returns the results in the structure
% FIT. The fits are performed using
% FITALLMODELSFOROBSERVATIONS2LAPLACE.  
%
% Optional arguments, provided as name-value pairs, are
%
% PNCBOTDATAFILE ['delayRegressData3']: The data file to use.
%
% Other optionals are passed to FITALLMODELSFOROBSERVATIONS2LAPLACE.

p = inputParser;
p.addOptional('pnCbotDataFile', 'delayRegressData3');
p.KeepUnmatched = true;
p.parse(varargin{:});

argsForFitFunction = UnpackStructureFieldsAsNameValuePairs(p.Unmatched);

Data = load(p.Results.pnCbotDataFile);

[mixtureInds, mixtureVals] = GetBinaryMixtureIndsForResponseLinearity();
whichMixture = mixtureVals(whichMixture,:);

indPure1 = GetIndexForBinaryMixtureConcentrationPair(whichMixture(1),0);
indPure2 = GetIndexForBinaryMixtureConcentrationPair(0,whichMixture(2));
indMix   = GetIndexForBinaryMixtureConcentrationPair(whichMixture(1), whichMixture(2));

Uall = squeeze(Data.pnCbot(whichCell,Data.whichBinsAll,[indPure1 indPure2 indMix],:));

offset = find(Data.whichBinsAll == Data.whichBinsToFit(1)) - 1;
whichBinsToFit = offset + (1:length(Data.whichBinsToFit));

Fit = FitAllModelsForObservations2Laplace(Uall(:,1:2,:), Uall(:,3,:),'whichBinsToFit',whichBinsToFit,argsForFitFunction{:});
Fit.tall = Data.binStarts(Data.whichBinsAll);
Fit.tfit = Data.binStarts(Data.whichBinsToFit);
Fit.whichCell = whichCell;
Fit.mixtureConcs = whichMixture;
Fit.indPure1 = indPure1;
Fit.indPure2 = indPure2;
Fit.indMix   = indMix;

