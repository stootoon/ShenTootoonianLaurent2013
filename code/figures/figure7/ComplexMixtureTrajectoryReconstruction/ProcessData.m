function ProcessData()
% function ProcessData()

global targetDir

figDir     = GetDataDirForFigure(7);
currDir    = GetCurrentDirFromPathString(fileparts(mfilename('fullpath')));
targetDir  = fullfile(figDir, currDir, 'recomputedData');

numMpCores = 6;

fprintf('Preparing datasets for fit...\n'); tic;
pnDataFile = PrepareDataForFits();
fprintf('Wrote "%s".\nDone in %1.1f seconds.\n', pnDataFile, toc); % ~ 3 secs

fprintf('Computing lasso weights...\n'); tic; 
lassoWeightsFile = ParComputeReconstructionWeightsUsingAdaptiveLasso(numMpCores, pnDataFile);
fprintf('Wrote "%s".\nDone in %1.1f seconds.\n', lassoWeightsFile, toc); % ~9000 seconds.

fprintf('Computing lasso weights for shuffled KCs...\n'); tic; 
lassoWeightsShuffledKcsFile = ParComputeLassoWeightsForShuffledKcs(numMpCores, pnDataFile, lassoWeightsFile);
fprintf('Wrote "%s".\nDone in %1.1f seconds.\n', lassoWeightsShuffledKcsFile, toc); % ~14000 seconds

fprintf('Computing lasso weights for shuffled PNs...\n'); tic; 
lassoWeightsShuffledPnsFile = ParComputeLassoWeightsForShuffledPns(numMpCores, pnDataFile); % ~190000 seconds
fprintf('Wrote "%s".\nDone in %1.1f seconds.\n', lassoWeightsShuffledPnsFile, toc); %

fprintf('Computing data for reconstructions...\n'); tic;
dataForPnRecFile = ComputeDataForPnReconstructions();
fprintf('Wrote "%s".\nDone in %1.1f seconds.\n', dataForPnRecFile, toc);

function outputFile = PrepareDataForFits()
% function outputFile = PrepareDataForFits()
%
% Prepares the dataset to be uesd for the fits, using binned PN and KC spike times.

global targetDir

pnSpt = LoadTocSpikeTimes('rawpn');
kcSpt = LoadTocSpikeTimes('rawkc');

t0 = 2.1;
t1 = 3.1;
binSize = 0.1;

Xpn = CountSpikesInBinsAndAverageAcrossTrials(pnSpt, {1:7}, 1:44, [1:174], 'startTime', t0, 'endTime', t1, 'binSize', binSize);
Ykc = CountSpikesInBinsAndAverageAcrossTrials(kcSpt, {1:7}, 1:44, [1:209], 'startTime', t0, 'endTime', t1, 'binSize', binSize);

[numPns,numBins,numOdors] = size(Xpn);
numKcs = size(Ykc,1);

% Reshape the data so that each column contains the binned responses of one cell.
% The rows are ordered as odor1 bin1...binN, odor2 bin1...binN, etc...
U = reshape(Xpn,[],numBins*numOdors)'; 
V = reshape(Ykc,[],numBins*numOdors)';

% Use only KCs that produced at least one responsive bin.
inz = find(sum(V));
V   = V(:,inz); 

numKcs     = size(V,2);
numSamples = numBins*numOdors;

fileName   = sprintf('data%d_t%1.1f_to_%1.1f', round(binSize*1000), t0, t1);
fileName   = [strrep(fileName, '.', '_') '.mat'];
outputFile = fullfile(targetDir, fileName);

save(outputFile, 'U', 'V');

function outputFile = ParComputeReconstructionWeightsUsingAdaptiveLasso(numMpCores, pnDataFile)
global targetDir

data = load(pnDataFile);
V = data.V;
U = data.U;

[numSamples, numKcs] = size(V);
[numSamples, numPns] = size(U);
whichPns = 1:numPns;
Results = [];

SetRandomSeed1(whichPns(1));

Results = zeros(numel(whichPns),numKcs+3,2); % +1 for offset weight, +2 for last stage, +3 for exitFlag
if (MatlabPoolWrapper('size')>0)
  MatlabPoolWrapper close;
end
MatlabPoolWrapper('open', numMpCores);

for i = 1:numel(whichPns)
  [b,b0,stats,exitFlag] = MultiStageAdaptiveLassoRegress(V,U(:,whichPns(i)),'numStages',10);
  Results(i,1:numKcs+1,1) = [b;b0];
  Results(i,1:numKcs+1,2) = [stats.B(:,1);stats.b0(1)];
  Results(i,numKcs+2,:) = stats.lastStage;
  Results(i,numKcs+3,:) = exitFlag;  
end

fileName = 'prunedWeightsAdaptiveLasso.mat';
outputFile = fullfile(targetDir, fileName);
save(outputFile, 'Results', 'pnDataFile');

function outputFile = ParComputeLassoWeightsForShuffledKcs(numMpCores, pnDataFile, lassoWeightsFile)
global targetDir

data = load(pnDataFile);
conn = load(lassoWeightsFile);

% The previously computed reconstruction weights
B = conn.Results(:,1:206,1)'; 
B0= conn.Results(:,207,1)';

U = data.U;
V = data.V;
[numSamples,numPns] = size(U);
[numSamples,numKcs] = size(V);

numOdors = length(GetOdorsList);
numBinsPerOdor = numSamples/numOdors;

bo_c2cbo = @(X) permute(reshape(X,numBinsPerOdor,numOdors,[]),[3 1 2]);
cbo2bo_c = @(Xcbo) reshape(Xcbo,[],numBinsPerOdor*numOdors)';

Vcbo = bo_c2cbo(V);
Ucbo = bo_c2cbo(U);

df = sum(B~=0,1); % The number of non-zero weights used in each reconstruction
Bsh = zeros(numKcs,numPns);
b0sh = zeros(1,numPns);

numKcShuffles = 100;
whichKcShuffles = 1:numKcShuffles;
Results = {};

SetRandomSeed1(whichKcShuffles(1));
if (MatlabPoolWrapper('size')>0)
  MatlabPoolWrapper close;
end
MatlabPoolWrapper('open', numMpCores);

for i = 1:numel(whichKcShuffles)
  % Shuffle the PNs
  VshCbo = CreateFakeCellsByShufflingOdorResponses(Vcbo);
  Vsh    = cbo2bo_c(VshCbo);
  Bsh    = zeros(numKcs,numPns);
  b0sh   = zeros(1,numPns);
  parfor j = 1:numPns
    [BB,bb,lambda,stats] = LassoRegress2(Vsh,U(:,j));
    dfsh = sum(BB~=0,1);
    % Find the reconstruction that used the same number of weights 
    indMatch = argmin(abs(dfsh-df(j))); 
    indMatch = indMatch(end); % lambda goes highest to lowest, so in case of a tie, pick the soln with the lower penalty
    Bsh(:,j) = BB(:,indMatch);
    b0sh(j)  = bb(indMatch);
  end
  Results{i} = {Vsh,Bsh,b0sh};
end

fileName   = 'lassoWeightsKcShuffle.mat'; 
outputFile = fullfile(targetDir, fileName);
save(outputFile, 'Results', 'pnDataFile', 'lassoWeightsFile');

function outputFile = ParComputeLassoWeightsForShuffledPns(numMpCores, pnDataFile)
global targetDir

data = load(pnDataFile);

V = data.V;
U = data.U;

[numSamples, numKcs] = size(V);
[numSamples, numPns] = size(U);

numOdors = length(GetOdorsList);
numBinsPerOdor = numSamples / numOdors;

bo_c2cbo = @(X) permute(reshape(X,numBinsPerOdor,numOdors,[]),[3 1 2]);
cbo2bo_c = @(X) reshape(X,[],numSamples)';

Vcbo = bo_c2cbo(V);
Ucbo = bo_c2cbo(U);

numPnShuffles   = 50;
whichPnShuffles = 1:numPnShuffles;

Results = {};

SetRandomSeed1(whichPnShuffles(1));
if (MatlabPoolWrapper('size')>0)
  MatlabPoolWrapper close;
end
MatlabPoolWrapper('open', numMpCores);

for i = 1:length(whichPnShuffles) 
  UshCbo = CreateFakeCellsByShufflingOdorResponses(Ucbo);
  Ush = cbo2bo_c(UshCbo);
  Bsh = zeros(numKcs,numPns);
  b0sh = zeros(1,numPns);
  exitFlags = zeros(1,numPns);
  for j = 1:numPns
    [Bsh(:,j),b0sh(j),stats,exitFlag(j)] = MultiStageAdaptiveLassoRegress(V,Ush(:,j),'numStages',20,'LassoFunction',@LassoRegress2);
  end
  Results{i} = {Ush, Bsh, b0sh, exitFlag};
end

fileName   = 'lassoWeightsPnShuffle.mat'; 
outputFile = fullfile(targetDir, fileName);
save(outputFile, 'Results', 'pnDataFile');


function outputFile = ComputeDataForPnReconstructions(sourceDir)
% Collects metrics for the trajectory reconstructions, specifically
% the correlation coefficients and the normalized SSE.
global targetDir

if (nargin==0)
  sourceDir = targetDir;
end

data     = load(fullfile(sourceDir, 'data100_t2_1_to_3_1.mat'));
conn     = load(fullfile(sourceDir, 'prunedWeightsAdaptiveLasso.mat'));
connPnSh = load(fullfile(sourceDir, 'lassoWeightsPnShuffle.mat'));
connKcSh = load(fullfile(sourceDir, 'lassoWeightsKcShuffle.mat'));

B = conn.Results(:,1:206,1)';
B0= conn.Results(:,207,1)';

U = data.U;
V = data.V;

[numSamples,numPns] = size(U);
[Numsamples,numKcs] = size(V);

numOdors = numel(GetOdorsList);

[UrecCbo,UpnCbo] = ReconstructPnTrajectories(V,U,B,B0);

[VkcSh, VkcShCbo, UkcShRec, UkcShRecCbo, BkcSh, b0kcSh] = CollectDataForShuffledKcs(connKcSh.Results);
[UpnSh, UpnShCbo, BpnSh, b0pnSh, UpnShRec, UpnShRecCbo] = CollectDataForShuffledPns(connPnSh.Results,V);

ssenData = {};
ccData   = {};
ssenData{1} = TrajectoryComparisonMetrics(UpnCbo, UrecCbo, 'ssen', @(M) mean(M));
ccData{1}   = TrajectoryComparisonMetrics(UpnCbo, UrecCbo, 'corr', @(M) mean(M));

numKcShuffles = size(UkcShRecCbo,4);
ssenData{2}   = zeros(numKcShuffles, numOdors);
ccData{2}     = zeros(numKcShuffles, numOdors);
for i = 1:numKcShuffles
  ssenData{2}(i,:) = TrajectoryComparisonMetrics(UpnCbo, UkcShRecCbo(:,:,:,i), 'ssen', @(M) mean(M));
  ccData{2}(i,:)   = TrajectoryComparisonMetrics(UpnCbo, UkcShRecCbo(:,:,:,i), 'corr', @(M) mean(M));
end
ssenData{2} = ssenData{2}(:);
ccData{2}   = ccData{2}(:);

numPnShuffles = size(UpnShRecCbo,4);
ssenData{3}   = zeros(numPnShuffles, numOdors);
ccData{3}     = zeros(numPnShuffles, numOdors);
for i = 1:numPnShuffles
  ssenData{3}(i,:) = TrajectoryComparisonMetrics(UpnShCbo(:,:,:,i), UpnShRecCbo(:,:,:,i), 'ssen', @(M) mean(M));
  ccData{3}(i,:)   = TrajectoryComparisonMetrics(UpnShCbo(:,:,:,i), UpnShRecCbo(:,:,:,i), 'corr', @(M) mean(M));
end
ssenData{3} = ssenData{3}(:);
ccData{3}   = ccData{3}(:);

medSsen = cellfun(@median, ssenData);
ubSsen  = cellfun(@(x) prctile(x, 95), ssenData);
lbSsen  = cellfun(@(x) prctile(x,  5), ssenData);

medCc = cellfun(@median, ccData);
ubCc  = cellfun(@(x) prctile(x, 95), ccData);
lbCc  = cellfun(@(x) prctile(x,  5), ccData);

fileName   = 'dataForPnRec.mat'; 
outputFile = fullfile(targetDir, fileName);
save(outputFile, 'ssenData', 'ccData', 'medSsen', 'ubSsen', 'lbSsen', 'medCc', 'ubCc', 'lbCc');
