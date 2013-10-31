function ProcessData()
% function ProcessData()

binaryMixturePairs = [0 140; 140 140; 140 0];

figDir    = GetDataDirForFigure(3);
currDir   = GetCurrentDirFromPathString(fileparts(mfilename('fullpath')));
targetDir = fullfile(figDir, currDir, 'recomputedData');

fileName = fullfile(targetDir, 'corrDistResults1to1Mixture.mat');

startTime = tic;
computeCorrelationDistancesForPairs(binaryMixturePairs, fileName);
fprintf('Wrote "%s" in %1.3f seconds.\n', fileName, toc(startTime));

concentrationSeriesPairs = [0   30;
       0   60;
       0   80;
       0  100;
       0  120;
       0  140;
       30  30;
       60  60;
       80  80;
       100 100;
       140 140
       30  0;
       60  0;
       80  0;
       100 0;
       120 0;
       140 0;];

fileName = fullfile(targetDir, 'corrDistResultsConcentrationSeries.mat');
startTime = tic;
computeCorrelationDistancesForPairs(concentrationSeriesPairs, fileName);
fprintf('Wrote "%s" in %1.3f seconds.\n', fileName, toc(startTime));

morphPairs = [0   140;
       30  140;
       60  140;
       80  140;
       100 140;
       120 140;
       140 140;
       140 120;
       140 100;
       140 80;
       140 60
       140 30;
       140 0]
fileName = 'corrDistResultsMixtureMorphs.mat';
startTime = tic;
computeCorrelationDistancesForPairs(morphPairs, fileName);
fprintf('Wrote "%s" in %1.3f seconds.\n', fileName, toc(startTime));


function computeCorrelationDistancesForPairs(whichPairs, fileName)
whichOdorInds = cell2mat(map(@(x) GetIndexForBinaryMixtureConcentrationPair(x(1),x(2)), whichPairs'))

X  = LoadTocSpikeTimes('rawpn_binary_mixtures');
Y1 = CountSpikesInBinsAndAverageAcrossTrials(X, 3:6,  whichOdorInds, 1:168, 'startTime', 2, 'endTime', 5, 'numAllTrials', 10, 'numAllOdors', 27);
Y2 = CountSpikesInBinsAndAverageAcrossTrials(X, 7:10, whichOdorInds, 1:168, 'startTime', 2, 'endTime', 5, 'numAllTrials', 10, 'numAllOdors', 27);

f  = @(M) min(M(:));
Y1 = reshape(Y1,168,[],numel(whichOdorInds));
Y2 = reshape(Y2,168,[],numel(whichOdorInds));
D  = ComputeCorrelationDistance(Y1,Y2);
Ds = SummarizeCorrelationDistanceMatrices(D,f);

save(fileName, 'D', 'Ds', 'whichPairs');
