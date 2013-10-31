function ProcessData()
% function ProcessData()

global targetDir
figDir    = GetDataDirForFigure(8);
currDir   = GetCurrentDirFromPathString(fileparts(mfilename('fullpath')));
targetDir = fullfile(figDir, currDir, 'recomputedData');

numMpCores = 6;
if (MatlabPoolWrapper('size')>0)
  MatlabPoolWrapper('close');
end
MatlabPoolWrapper('open', numMpCores);

disp('Preparing data for classification...');
startTime = tic;
dataFile = PrepareData(0, 4.5, ''); % ~60 seconds
fprintf('Wrote "%s" in %1.1f secs.\n', dataFile, toc(startTime));

disp('Preparing data for classification (short)...');
startTime = tic;
dataFileShort = PrepareData(1.9, 4.4, 'Short'); % ~60 seconds
fprintf('Wrote "%s" in %1.1f secs.\n', dataFileShort, toc(startTime));

disp('Generating PN shuffle Ids for categorization and generalization...');
startTime = tic;
pnShuffleIdsFile = GeneratePnShuffleIds();
fprintf('Wrote "%s" in %1.1f secs.\n', pnShuffleIdsFile, toc(startTime));

%% ***** Identity Decoding Performance **********
disp('Computing identity decoding results...');
startTime = tic;
identityDecodingFile = ComputeDataForIdentityDecoding();
fprintf('Wrote "%s" in %1.1f secs.\n', identityDecodingFile, toc(startTime));

disp('Computing identity decoding results for top PNs...');
aucResultsFile = fullfile(GetDataDirForFigure(6),'PnKcAuc_noduplicates','originalData','Mfig_noduplicates.mat');
startTime = tic;
identityDecodingTopPnsFile = ComputeDataForIdentityDecodingTopPns(aucResultsFile); % 103.0 secs
fprintf('Wrote "%s" in %1.1f secs.\n', identityDecodingTopPnsFile, toc(startTime));

disp('Computing identity for random PNs...');
startTime = tic;
identityDecodingRandomPnsFile = ParComputeIdentityDecodingResultsForRandomPns(dataFile); % ~10268 secs
fprintf('Wrote "%s" in %1.1f secs.\n', identityDecodingRandomPnsFile, toc(startTime));

%% ***** Categorization Performance *************
disp('Computing categorization performance for PNs...');
startTime = tic;
pnCatResultsFile = ParComputeCategorizationResultsForPns(dataFileShort); % ~ 45 mins
fprintf('Wrote "%s" in %1.1f secs.\n', pnCatResultsFile, toc(startTime));

disp('Computing categorization performance for KCs...');
startTime = tic;
kcCatResultsFile = ParComputeCategorizationResultsForKcs(dataFileShort); % ~45 mins
fprintf('Wrote "%s" in %1.1f secs.\n', kcCatResultsFile, toc(startTime));

disp('Computing categorization performance for top PNs...');
startTime = tic;
pnCatResultsTopPnsFile = ParComputeCategorizationResultsForTopPns(dataFileShort, identityDecodingTopPnsFile); % ~45 mins
fprintf('Wrote "%s" in %1.1f secs.\n', pnCatResultsTopPnsFile, toc(startTime));

disp('Computing categorization peformance for random PNs...');
startTime = tic;
pnCatResultsRandomPnsFile = ParComputeCategorizationResultsForRandomPns(dataFileShort, pnShuffleIdsFile); % ~220000 seconds for 100 shuffles
fprintf('Wrote "%s" in %1.1f secs.\n', pnCatResultsRandomPnsFile, toc(startTime));

%% ***** Generalization Performance *************
disp('Computing generalization performance for PNs...');
startTime = tic;
pnGenResultsFile = ParComputeGeneralizationResultsForPns(dataFileShort); % ~541 secs 
fprintf('Wrote "%s" in %1.1f secs.\n', pnGenResultsFile, toc(startTime));

disp('Computing generalization performance for KCs...');
startTime = tic;
kcGenResultsFile = ParComputeGeneralizationResultsForKcs(dataFileShort); % ~605 secs
fprintf('Wrote "%s" in %1.1f secs.\n', kcGenResultsFile, toc(startTime));

disp('Computing generalization performance for top PNs...');
startTime = tic;
pnGenResultsTopPnsFile = ParComputeGeneralizationResultsForTopPns(dataFileShort, identityDecodingTopPnsFile); % ~400 secs
fprintf('Wrote "%s" in %1.1f secs.\n', pnGenResultsTopPnsFile, toc(startTime));

disp('Computing generalization peformance for random PNs...');
startTime = tic;
pnGenResultsRandomPnsFile = ParComputeGeneralizationResultsForRandomPns(dataFileShort, pnShuffleIdsFile); % ~34000 seconds
fprintf('Wrote "%s" in %1.1f secs.\n', pnGenResultsRandomPnsFile, toc(startTime));


function outputFileName = PrepareData(t0, t1, suffix)
global targetDir

pnSpt = ConvertSpikeTimesFromSparseToFull(LoadTocSpikeTimes('rawpn'));
kcSpt = ConvertSpikeTimesFromSparseToFull(LoadTocSpikeTimes('rawkc'));

startTime = t0;
endTime   = t1;
binSize   = 0.1;
stepSize  = 0.025;
overlap   = binSize - stepSize;

odors     = GetOdorsList;
numTrials = 7;
numOdors  = length(odors);

numPns = size(pnSpt,2)/numTrials/numOdors;
numKcs = size(kcSpt,2)/numTrials/numOdors;

disp('Loading PN data...');
[PnCbot, pnBinStarts, pnImove] = CountSpikesInBinsAndAverageAcrossTrials2(pnSpt, {1,2,3,4,5,6,7}, 1:numOdors, 1:numPns, 'startTime',startTime,'endTime',endTime,'binSize',binSize,'overlap',overlap);
fprintf('  %d/%d spikes were moved to accommodate bin edges.\n', nnz(pnImove), sum(~isinf(pnSpt(:))));

disp('Loading KC data...');
[KcCbot, kcBinStarts, kcImove] = CountSpikesInBinsAndAverageAcrossTrials2(kcSpt, {1,2,3,4,5,6,7}, 1:numOdors, 1:numKcs, 'startTime',startTime,'endTime',endTime,'binSize',binSize,'overlap',overlap);
fprintf('  %d/%d spikes were moved to accommodate bin edges.\n', nnz(kcImove), sum(~isinf(kcSpt(:))));

outputFileName = fullfile(targetDir, ['dataForClassification' suffix]);
save(outputFileName, 'PnCbot', 'KcCbot', 'pnBinStarts', 'kcBinStarts', 'startTime', 'endTime', 'binSize', 'stepSize', 'overlap', 'pnImove', 'kcImove');

function outputFile = ComputeDataForIdentityDecoding()
global targetDir
%% Compute the stats for the raw data
lambda = 1;
Data = load(fullfile(targetDir, 'dataForClassification.mat'));
[ResultsPn, ScoresPn] = ClassifyOdorIdentityUsingRlsc(Data.PnCbot, lambda);
[ResultsKc, ScoresKc] = ClassifyOdorIdentityUsingRlsc(Data.KcCbot, lambda);

outputFile = fullfile(targetDir, 'dataIdentityDecoding.mat');
save(outputFile, 'Data', 'lambda', 'ResultsPn', 'ResultsKc', 'ScoresPn', 'ScoresKc');

function outputFile = ComputeDataForIdentityDecodingTopPns(aucResultsFile)
global targetDir
%% Compute the stats using the top PNs for each odor
lambda = 1;
% Load the AUC data
AucData = load(aucResultsFile);
AucData = AucData.M;
% Find the top PN for each odor.

cells = [AucData.RocPn.cell];
cmps  = [AucData.RocPn.cmp];
auc   = [AucData.RocPn.auc];

topPnStructIds = arrayfun(@(i) argmax((cmps==i).*auc,1),1:8);
topPnIds = arrayfun(@(i) cells(argmax((cmps==i).*auc,1)),1:8);
topPnAuc = auc(topPnStructIds);
Data  = load(fullfile(targetDir, 'dataForClassification.mat'));
Xcbot = Data.PnCbot;
Xcbot = Xcbot(topPnIds,:,:,:);
[ResultsTopPns, ScoresTopPns] = ClassifyOdorIdentityUsingRlsc(Xcbot, lambda);

outputFile = fullfile(targetDir, 'dataIdentityDecodingTopPns.mat');
save(outputFile, 'topPnIds', 'topPnAuc', 'ResultsTopPns', 'ScoresTopPns'); 

function outputFile = ParComputeIdentityDecodingResultsForRandomPns(dataFile)
global targetDir
lambda           = 1; % Regularization parameter
numPnsPerShuffle = 8;
numShuffles      = 100;
whichShuffles    = 1:numShuffles;
Results          = {};

rndSeed = SetRandomSeed1(whichShuffles(1),'override','yes');
Data    = load(dataFile);
PnCbot  = Data.PnCbot;
numPns  = size(PnCbot,1);

for i = 1:numel(whichShuffles)
  whichPns         = first(randperm(numPns), numPnsPerShuffle);
  Xcbot            = Data.PnCbot(whichPns,:,:,:);
  [Output, Scores] = ClassifyOdorIdentityUsingRlsc(Xcbot, lambda);
  Results{i}       = {Output, Scores, whichPns, rndSeed};
end

outputFile = fullfile(targetDir, 'dataIdentityDecodingPnRandomShuffles.mat');
save(outputFile, 'Results');

function outputFile = ParComputeCategorizationResultsForPns(dataFileShort)
global targetDir
lambda = 1;

cmps = 'ABCDWXYZ';
numRuns = length(cmps);
whichRuns = 1:numRuns;
numReps = 20;
verbose = 0;
Results = {};

SetRandomSeed1(whichRuns(1));

Data = load(dataFileShort);

addpath(fullfile(GetRootDir,'code','figures','figure8'));

for i = 1:numel(whichRuns)
  whichCmp = whichRuns(i);
  [Output, OdorIds] = ComputeCategorizationPerformanceForComponentFast(Data.PnCbot, cmps(whichCmp), numReps, lambda, verbose);
  Results{i} = {Output, OdorIds, cmps(whichCmp)};
end

outputFile = fullfile(targetDir, 'dataCategoryDecodingPns.mat');
save(outputFile, 'Results');

function outputFile = ParComputeCategorizationResultsForKcs(dataFileShort)
global targetDir

lambda = 1;

cmps = 'ABCDWXYZ';
numRuns = length(cmps);
whichRuns = 1:numRuns;
numReps = 20;
verbose = 0;
Results = {};

SetRandomSeed1(whichRuns(1));

Data = load(dataFileShort);

addpath(fullfile(GetRootDir,'code','figures','figure8'));

for i = 1:numel(whichRuns)
  whichCmp = whichRuns(i);
  [Output, OdorIds] = ComputeCategorizationPerformanceForComponentFast(Data.KcCbot, cmps(whichCmp), numReps, lambda, verbose);
  Results{i} = {Output, OdorIds, cmps(whichCmp)};
end

outputFile = fullfile(targetDir, 'dataCategoryDecodingKcs.mat');
save(outputFile, 'Results');

function outputFile = GeneratePnShuffleIds()
global targetDir

SetRandomSeed1(1,'override','yes');

numPns           = 174;
numPnsPerShuffle = 8;
numShuffles      = 1000;
[foo, IpnSh]     = sort(rand(174,numShuffles));
IpnSh            = IpnSh(1:numPnsPerShuffle,:);

outputFile = fullfile(targetDir, 'pnShuffleIds.mat');
save(outputFile, 'IpnSh');

function outputFile = ParComputeCategorizationResultsForTopPns(dataFileShort, topPnsFile)
global targetDir

addpath(fullfile(GetRootDir,'code','figures','figure8'));

lambda = 1; % Regularization parameter
DataTop = load(topPnsFile);
cmps      = 'ABCDWXYZ';
numRuns   = length(cmps);
whichRuns = 1:numRuns;
numReps = 20;
verbose = 0;
Results = {};

SetRandomSeed1(whichRuns(1));

Data = load(dataFileShort);
for i = 1:numel(whichRuns)
  whichCmp = whichRuns(i);
  Xcbot = Data.PnCbot(DataTop.topPnIds,:,:,:);
  Output = ComputeCategorizationPerformanceForComponentFast(Xcbot, cmps(whichCmp), numReps, lambda, verbose);
  Results{i} = {Output, cmps(whichCmp)};
end

outputFile = fullfile(targetDir, 'dataCategoryDecodingTopPns.mat');
save(outputFile, 'Results');

function outputFile = ParComputeCategorizationResultsForRandomPns(dataFileShort, pnShuffleIdsFile)
global targetDir

addpath(fullfile(GetRootDir,'code','figures','figure8'));

lambda      = 1; % Regularization parameter
numShuffles = 10; % A smaller number because it takes a long time to do the full 100.
cmps        = 'ABCDWXYZ';
numRuns = numShuffles*length(cmps);
whichRuns = 1:numRuns;
numReps = 20;
verbose = 0;
Results = {};

SetRandomSeed1(whichRuns(1));
Data       = load(dataFileShort);
PnCbot     = Data.PnCbot;
ShuffleIds = load(pnShuffleIdsFile);
for i = 1:numel(whichRuns)
  whichShuffle = floor((whichRuns(i)-1)/length(cmps))+1;
  whichCmp = mod(whichRuns(i)-1,length(cmps))+1;
  whichPns = ShuffleIds.IpnSh(:,whichShuffle);
  Xcbot = Data.PnCbot(whichPns,:,:,:);
  Output = ComputeCategorizationPerformanceForComponentFast(Xcbot, cmps(whichCmp), numReps, lambda, verbose);
  Results{i} = {Output, whichShuffle, cmps(whichCmp), whichPns};
end

outputFile = fullfile(targetDir, 'dataCategoryDecodingPnRandomShuffles.mat');
save(outputFile, 'Results','-v7.3');

function outputFile = ParComputeGeneralizationResultsForPns(dataFileShort)
global targetDir

addpath(fullfile(GetRootDir,'code','figures','figure8'));

lambda = 1; % Regularization parameter

cmps    = 'ABCDWXYZ';
numRuns = length(cmps);
whichRuns = 1:numRuns;
numReps = 20;
verbose = 0;
Results = {};

SetRandomSeed1(whichRuns(1));

Data = load(dataFileShort);

for i = 1:numel(whichRuns)
  whichCmp = whichRuns(i);
  Output = ComputeGeneralizationPerformanceForComponent(Data.PnCbot, cmps(whichCmp), numReps, lambda, verbose);
  Results{i} = {Output, cmps(whichCmp)};
end

outputFile = fullfile(targetDir, 'dataGeneralizationDecodingPns.mat');
save(outputFile, 'Results');

function outputFile = ParComputeGeneralizationResultsForKcs(dataFileShort)
global targetDir

addpath(fullfile(GetRootDir,'code','figures','figure8'));

lambda = 1; % Regularization parameter

cmps    = 'ABCDWXYZ';
numRuns = length(cmps);
whichRuns = 1:numRuns;
numReps = 20;
verbose = 0;
Results = {};

SetRandomSeed1(whichRuns(1));

Data = load(dataFileShort);

for i = 1:numel(whichRuns)
  whichCmp = whichRuns(i);
  Output = ComputeGeneralizationPerformanceForComponent(Data.KcCbot, cmps(whichCmp), numReps, lambda, verbose);
  Results{i} = {Output, cmps(whichCmp)};
end

outputFile = fullfile(targetDir, 'dataGeneralizationDecodingKcs.mat');
save(outputFile, 'Results');

function outputFile = ParComputeGeneralizationResultsForTopPns(dataFileShort, topPnsFile)
global targetDir

addpath(fullfile(GetRootDir,'code','figures','figure8'));

lambda = 1; % Regularization parameter

DataTop = load(topPnsFile);

cmps = 'ABCDWXYZ';
numRuns = length(cmps);
whichRuns = 1:numRuns;
numReps = 20;
verbose = 0;
Results = {};

SetRandomSeed1(whichRuns(1));

Data = load(dataFileShort);

for i = 1:numel(whichRuns)
  whichCmp = i;
  Xcbot = Data.PnCbot(DataTop.topPnIds,:,:,:);
  Output = ComputeGeneralizationPerformanceForComponent(Xcbot, cmps(whichCmp), numReps, lambda, verbose);
  Results{i} = {Output, cmps(whichCmp)};
end

outputFile = fullfile(targetDir, 'dataGeneralizationDecodingTopPns.mat');
save(outputFile, 'Results');

function outputFile = ParComputeGeneralizationResultsForRandomPns(dataFileShort, pnShuffleIdsFile)
global targetDir

addpath(fullfile(GetRootDir,'code','figures','figure8'));

lambda = 1; % Regularization parameter
numShuffles = 100;
cmps = 'ABCDWXYZ';
numRuns = numShuffles*length(cmps);
whichRuns = 1:numRuns;
numReps = 20;
verbose = 0;
Results = {};

SetRandomSeed1(whichRuns(1));

Data = load(dataFileShort);

PnCbot = Data.PnCbot;

ShuffleIds = load(pnShuffleIdsFile);

for i = 1:numel(whichRuns)
  whichShuffle = floor((whichRuns(i)-1)/length(cmps))+1;
  whichCmp = mod(whichRuns(i)-1,length(cmps))+1;
  whichPns = ShuffleIds.IpnSh(:,whichShuffle);
  Xcbot = Data.PnCbot(whichPns,:,:,:);
  Output = ComputeGeneralizationPerformanceForComponent(Xcbot, cmps(whichCmp), numReps, lambda, verbose);
  Results{i} = {Output, whichShuffle, cmps(whichCmp), whichPns};
end

outputFile = fullfile(targetDir, 'dataGeneralizationDecodingPnRandomShuffles.mat');
save(outputFile, 'Results');
