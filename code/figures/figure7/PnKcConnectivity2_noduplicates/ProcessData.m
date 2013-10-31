function ProcessData
global targetDir
figDir     = GetDataDirForFigure(7);
currDir    = GetCurrentDirFromPathString(fileparts(mfilename('fullpath')));
targetDir  = fullfile(figDir, currDir, 'recomputedData');

numMpCores = 6;

fprintf('Preparing datasets for fit...\n'); tic;
pnDataFile = PrepareDataForFits();
fprintf('Wrote "%s".\nDone in %1.1f seconds.\n', pnDataFile, toc); % ~ 6 secs

fprintf('Computing lasso weights...\n'); tic; 
lassoWeightsFile = ParComputePnToKcConnUsingAdaptiveLasso(numMpCores, pnDataFile); % ~270 seconds
fprintf('Wrote "%s".\nDone in %1.1f seconds.\n', lassoWeightsFile, toc); % ~9000 seconds.

fprintf('Computing lasso weights for shuffled KCs...\n'); tic; 
lassoWeightsShuffledKcsFile = ParComputeLassoWeightsForShuffledKcs(numMpCores, pnDataFile, lassoWeightsFile);
fprintf('Wrote "%s".\nDone in %1.1f seconds.\n', lassoWeightsShuffledKcsFile, toc); % ~43000 seconds

fprintf('Computing lasso weights for shuffled PNs...\n'); tic; 
lassoWeightsShuffledPnsFile = ParComputeLassoWeightsForShuffledPns(numMpCores, pnDataFile, lassoWeightsFile); % ~6540 seconds
fprintf('Wrote "%s".\nDone in %1.1f seconds.\n', lassoWeightsShuffledPnsFile, toc); %

fprintf('Computing data for reconstructions...\n'); tic;
dataForKcRecFile = ComputeDataForKcReconstructions(pnDataFile, lassoWeightsFile, lassoWeightsShuffledPnsFile, lassoWeightsShuffledKcsFile); % ~14 seconds
fprintf('Wrote "%s".\nDone in %1.1f seconds.\n', dataForKcRecFile, toc);

function outputFile = PrepareDataForFits()
global targetDir

binSize = 0.1;
t0      = 2.1;
t1      = 3.1;

pnSpt = ConvertSpikeTimesFromSparseToFull(LoadTocSpikeTimes('rawpn'));
U     = CountSpikesInBinsAndAverageAcrossTrials(pnSpt,1:7,1:44,1:174,'startTime',t0,'endTime',t1,'binSize',binSize);
U     = reshape(U,size(U,1),[])';

kcSpt = ConvertSpikeTimesFromSparseToFull(LoadTocSpikeTimes('rawkc'));
V     = CountSpikesInBinsAndAverageAcrossTrials(kcSpt,1:7,1:44,1:209,'startTime',t0,'endTime',t1,'binSize',binSize);
V     = reshape(V,size(V,1),[])';
V     = V(:,find(sum(V))); % Grab the KCs that have at least one responding bin.

fileName   = sprintf('data%d_t%1.1f_to_%1.1f', round(binSize*1000), t0, t1);
fileName   = [strrep(fileName, '.', '_') '.mat'];
outputFile = fullfile(targetDir, fileName);
save(outputFile,'U','V');

function outputFile = ParComputePnToKcConnUsingAdaptiveLasso(numMpCores, pnDataFile)
global targetDir

data = load(pnDataFile);
V = data.V;
U = data.U;

[numSamples, numKcs] = size(V);
[numSamples, numPns] = size(U);

whichKcs = 1:numKcs;
Results  = [];

LASSOF   = @LassoRegress2NonNegative;

SetRandomSeed1(whichKcs(1));

Results = zeros(numel(whichKcs),numPns+3,2); % +1 for offset weight, +2 for last stage, +3 for exitFlag

if (MatlabPoolWrapper('size')>0)
  MatlabPoolWrapper close;
end
MatlabPoolWrapper('open', numMpCores);

for i = 1:numel(whichKcs)
  [b,b0,stats,exitFlag]   = MultiStageAdaptiveLassoRegress(U,V(:,whichKcs(i)),'numStages',10,'LassoFunction',LASSOF);
  Results(i,1:numPns+1,1) = [b;b0];
  Results(i,1:numPns+1,2) = [stats.B(:,1);stats.b0(1)];
  Results(i,numPns+2,:)   = stats.lastStage;
  Results(i,numPns+3,:)   = exitFlag;  
end

fileName = 'prunedWeightsAdaptiveLasso.mat';
outputFile = fullfile(targetDir, fileName);
save(outputFile, 'Results', 'pnDataFile');

function outputFile = ParComputeLassoWeightsForShuffledKcs(numMpCores, pnDataFile, lassoWeightsFile)
global targetDir

data    = load(pnDataFile);
weights = load(lassoWeightsFile);

U = data.U;
V = data.V;

[numSamples, numPns] = size(U);
[numSamples, numKcs] = size(V);

numOdors       = length(GetOdorsList);
numBinsPerOdor = numSamples/numOdors;

bo_c2cbo = @(X)    permute(reshape(X,numBinsPerOdor,numOdors,[]),[3 1 2]);
cbo2bo_c = @(Xcbo) reshape(Xcbo,[],numSamples)';

Ucbo = bo_c2cbo(U);
Vcbo = bo_c2cbo(V);

B = weights.Results(:,1:174,1)';
b0= weights.Results(:,175,1)';
df = sum(B~=0,1);

numShuffles = 250;
whichShuffles = 1:numShuffles;
Results = {};

SetRandomSeed1(whichShuffles(1));
if (MatlabPoolWrapper('size')>0)
  MatlabPoolWrapper close;
end
MatlabPoolWrapper('open', numMpCores);

for i = 1:length(whichShuffles)
  VshCbo = CreateFakeCellsByShufflingOdorResponses(Vcbo);
  Vsh    = cbo2bo_c(VshCbo);
  Bsh    = 0*B;
  b0sh   = 0*b0;
  for j = 1:numKcs % Can't use PARFOR - it's used inside MULTISTAGE...
    [Bsh(:,j),b0sh(j),stats,exitFlag] = MultiStageAdaptiveLassoRegress(U,Vsh(:,j),'numStages',1,'LassoFunction',@LassoRegress2NonNegative); % Do one stage - just the lasso.
  end
  Results{i} = {Vsh, Bsh, b0sh, stats, exitFlag};
end

fileName   = 'lassoWeightsKcShuffle.mat'; 
outputFile = fullfile(targetDir, fileName);
save(outputFile, 'Results', 'pnDataFile', 'lassoWeightsFile');

function outputFile = ParComputeLassoWeightsForShuffledPns(numMpCores, pnDataFile, lassoWeightsFile)
global targetDir

data    = load(pnDataFile);
weights = load(lassoWeightsFile);

U = data.U;
V = data.V;

[numSamples, numPns] = size(U);
[numSamples, numKcs] = size(V);

numOdors = length(GetOdorsList);
numBinsPerOdor = numSamples/numOdors;

bo_c2cbo = @(X) permute(reshape(X,numBinsPerOdor,numOdors,[]),[3 1 2]);
cbo2bo_c = @(Xcbo) reshape(Xcbo,[],numSamples)';

Ucbo = bo_c2cbo(U);
Vcbo = bo_c2cbo(V);

B  = weights.Results(:,1:174,1)';
b0 = weights.Results(:,175,1)';
df = sum(B~=0,1);

numShuffles = 1000;
whichShuffles = 1:numShuffles;
Results = {};

SetRandomSeed1(whichShuffles(1));
if (MatlabPoolWrapper('size')>0)
  MatlabPoolWrapper close;
end
MatlabPoolWrapper('open', numMpCores);

for i = 1:length(whichShuffles)
  UshCbo = CreateFakeCellsByShufflingOdorResponses(Ucbo);
  Ush = cbo2bo_c(UshCbo);
  Bsh = 0*B;
  b0sh = 0*b0;
  parfor j = 1:numKcs
    [W,w0,lambda,stats] = LassoRegress2NonNegative(Ush,V(:,j));
    ind = argmin(abs(stats.df-df(j)));
    if (length(ind)>1)
      ind = ind(argmin(lambda(ind))); % Take the one with a lower penalty.
    end
    Bsh(:,j) = W(:,ind);
    b0sh(j)  = w0(ind);
  end
  Results{i} = {Ush, Bsh, b0sh};
end

fileName   = 'lassoWeightsPnShuffle.mat'; 
outputFile = fullfile(targetDir, fileName);
save(outputFile, 'Results', 'pnDataFile', 'lassoWeightsFile');

function outputFile = ComputeDataForKcReconstructions(pnDataFile, lassoWeightsFile, lassoWeightsPnShuffleFile, lassoWeightsKcShuffleFile)
global targetDir

data        = load(pnDataFile);
weights     = load(lassoWeightsFile);
weightsPnSh = load(lassoWeightsPnShuffleFile);
weightsKcSh = load(lassoWeightsKcShuffleFile);

Results = weights.Results;
[UshPn, BshPn, b0ShPn] = CollectPnShufflingResults(weightsPnSh.Results);
[VshKc, BshKc, b0ShKc] = CollectKcShufflingResults(weightsKcSh.Results);

U = data.U;
V = data.V;

[numSamples, numPns] = size(U);
[numSamples, numKcs] = size(V);

B = weights.Results(:,1:174,1)';
b0= weights.Results(:,175,1)';

muB   = mean(B>0);
muMuB = mean(muB(muB>0))*size(B,1);
sdMuB = std(muB(muB>0));
seMuB = sdMuB/sqrt(sum(muB>0))*size(B,1);
maxMuB = max(muB)*size(B,1);

fprintf('Fraction not fit = %1.3f. Mean+/-SEM connectivity of remainder: %1.3f, %1.3f. Max: %1.3f\n',mean(muB==0), muMuB, seMuB, maxMuB);

Vrec = bsxfun(@plus, U*B, b0);
ssenRec = KcReconstructionComparisonMetrics(V,Vrec,'ssen',@Identity);
ccRec  = KcReconstructionComparisonMetrics(V,Vrec,'corr',@Identity);

numPnShuffles = size(UshPn,3);
numKcShuffles = size(VshKc,3);

VrecPnSh = zeros(numSamples,numKcs,numPnShuffles);
ssenRecPnSh = zeros(numPnShuffles,numKcs);
ccRecPnSh  = zeros(numPnShuffles,numKcs);

for i = 1:numPnShuffles
  % Reconstruct the real KCs using the shuffled PN data.
  VrecPnSh(:,:,i) = bsxfun(@plus, UshPn(:,:,i)*BshPn(:,:,i), b0ShPn(i,:));
  ssenRecPnSh(i,:) = KcReconstructionComparisonMetrics(V,VrecPnSh(:,:,i),'ssen',@Identity);
  ccRecPnSh(i,:) = KcReconstructionComparisonMetrics(V,VrecPnSh(:,:,i),'corr',@Identity);
end

VrecKcSh = zeros(numSamples,numKcs,numPnShuffles);
ssenRecKcSh = zeros(numKcShuffles,numKcs);
ccRecKcSh  = zeros(numKcShuffles,numKcs);
for i = 1:numKcShuffles
  % Reconstruct the shuffled KCs using real PN data.
  VrecKcSh(:,:,i) = bsxfun(@plus, U*BshKc(:,:,i), b0ShKc(i,:));
  ssenRecKcSh(i,:) = KcReconstructionComparisonMetrics(VshKc(:,:,i),VrecKcSh(:,:,i),'ssen',@Identity);
  ccRecKcSh(i,:) = KcReconstructionComparisonMetrics(VshKc(:,:,i),VrecKcSh(:,:,i),'corr',@Identity);  
end

% Define a KC to be well reconstructed if the SSEn of its
% reconstruction is less than the minimum if the PNs are shuffled, and
% the CC of its reconstruction is greater than the maximum if the PNs
% are shuffled.
indSig = find((ssenRec < min(ssenRecPnSh)) & (ccRec>max(ccRecPnSh)));

fileName   = 'dataForKcRec.mat'; 
outputFile = fullfile(targetDir, fileName);
save(outputFile, 'indSig', 'ssenRec', 'ssenRecPnSh', 'ssenRecKcSh', 'ccRec', 'ccRecPnSh', 'ccRecKcSh');
