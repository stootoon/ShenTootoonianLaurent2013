function ProcessData()
% function ProcessData()
global targetDir;

figDir    = GetDataDirForFigure(2);
thisDir   = GetCurrentDirFromPathString(fileparts(mfilename('fullpath')));
targetDir = fullfile(figDir, thisDir, 'recomputedData')

numMpCores = 8;

disp('Setting up data for regression...'); startTime = tic; % 22 seconds
dataFileForRegression = PrepareDataForRegression();
fprintf('Wrote "%s" in %1.3f seconds.\n\n', dataFileForRegression, toc(startTime));

disp('Computing model fits...'); startTime = tic; % 351 seconds
modelFitResultsFile   = ParComputeBestModels3Bayes1Laplace(numMpCores, dataFileForRegression);
fprintf('Wrote "%s" in %1.3f seconds.\n\n', modelFitResultsFile, toc(startTime));

disp('Computing colormaps...'); startTime = tic; % 84 seconds
colorMapResultsFile   = ComputeColormaps(modelFitResultsFile);
fprintf('Wrote "%s" in %1.3f seconds.\n\n', colorMapResultsFile, toc(startTime));

disp('Computing SNR Angle data...'); startTime = tic;
snrAngleFile = ComputeSnrAngleData(modelFitResultsFile);
fprintf('Wrote "%s" in %1.3f seconds.\n\n', snrAngleFile, toc(startTime));

function outputFileName = PrepareDataForRegression()
global targetDir
% PrepareDataForRegression()
t0      = 0;
t1      = 6;
t1end   = 13; % We need this when computing the SNRs, to get additional baseline data
binSize = 0.05;
pnSpt   = LoadTocSpikeTimes('rawpn_binary_mixtures');
pnCbot  = CountSpikesInBinsAndAverageAcrossTrials(pnSpt,arrayfun(@Identity, 1:10,'UniformOutput',false),1:27,1:168,'startTime',t0,'endTime',t1end,'binSize',binSize,'numAllTrials',10,'numAllOdors',27);

pureInd1   = GetIndexForBinaryMixtureConcentrationPair(60, 0);
pureInd2   = GetIndexForBinaryMixtureConcentrationPair(0,  140);
mixtureInd = GetIndexForBinaryMixtureConcentrationPair(60, 140);
cellInd    = 18;

numBins        = size(pnCbot,2);
binStarts      = (0:numBins-1)*binSize + t0;
t0fit          = 1.5;
t1fit          = 5;
whichBinsToFit = find(binStarts>=t0fit & binStarts<t1fit);
whichBinsAll   = find(binStarts>=t0    & binStarts<t1);

outputFileName = fullfile(targetDir, 'delayRegressData.mat');
save(outputFileName, 'pnCbot', 'binSize', 't0', 't1', 't1end', 'pureInd1', 'pureInd2', 'mixtureInd', 'cellInd', 'numBins', 'binStarts', 't0fit', 't1fit', 'whichBinsToFit', 'whichBinsAll');

function outputFileName = ParComputeBestModels3Bayes1Laplace(numMpCores, dataFile)
% Same as PARCOMPUTEBESTMODELS2, but uses Bayesian model selection,
% with gaussian prior on regression coefficients and lags.
%
% dataFile = 'delayRegressData.mat';
global targetDir

Data = load(dataFile);

pnCbot = Data.pnCbot;

numTrials      = size(pnCbot,4);
whichBinsAll   = Data.whichBinsAll;
whichBinsToFit = Data.whichBinsToFit;

offset = find(whichBinsAll == whichBinsToFit(1)) - 1;

[mixtureInds, mixtureVals] = GetBinaryMixtureIndsForResponseLinearity;

[numCells, numTotalBins, numOdors, numTrials] = size(pnCbot);
numMixtures = size(mixtureVals,1);

whichCells    = 1:numCells;
whichMixtures = 1:numMixtures;
params        = CartesianProduct({whichCells, whichMixtures});
numParams     = size(params,1);
whichParams   = 1:numParams;

% Model priors
logModelPriors = @(whichInputs) -sum(whichInputs)*log(numel(whichInputs));
sigmaLag = 1;
sigmaReg = 1;
lagLimit = 3;
variancePriorAlpha = 1;
variancePriorBeta  = 1;
fitOptions = {'whichBinsToFit',offset + (1:length(whichBinsToFit)),'sigmaLag',sigmaLag,'sigmaReg',sigmaReg, 'lagLimit', lagLimit, 'logModelPriors', logModelPriors, 'variancePriorAlpha', variancePriorAlpha, 'variancePriorBeta', variancePriorBeta};

Results = {};

% Start the computation
if (MatlabPoolWrapper('size')>0)
  MatlabPoolWrapper close
end

MatlabPoolWrapper('open', numMpCores);

numParamsToDo = numel(whichParams);

warning off;

Results = cell(numParamsToDo,1);

parfor i = 1:numParamsToDo
  thisParam    = params(whichParams(i),:);  
  whichCell    = thisParam(1);
  whichMixture = thisParam(2);
  
  indPure1 = GetIndexForBinaryMixtureConcentrationPair(mixtureVals(whichMixture,1),0);
  indPure2 = GetIndexForBinaryMixtureConcentrationPair(0,mixtureVals(whichMixture,2));
  indMix   = GetIndexForBinaryMixtureConcentrationPair(mixtureVals(whichMixture,1), mixtureVals(whichMixture,2));
    
  % Get the data 
  Xobs = squeeze(pnCbot(whichCell,whichBinsAll,[indPure1 indPure2],:));
  yobs = squeeze(pnCbot(whichCell,whichBinsAll,indMix,:));
  
  F = FitAllModelsForObservations2Laplace(Xobs, yobs, fitOptions{:});
  
  M = F.Models{F.bestModel};
  M2= F.Models{15};
  
  Results{i,1} = struct;
  
  Results{i,1}.bestModel = F.bestModel;
  Results{i,1}.descr = M.descr;
  Results{i,1}.indInputs = M.indInputs;
  Results{i,1}.w = M.w;
  Results{i,1}.v = M.v;
  Results{i,1}.lags = M.lags;
  
  Results{i,1}.r2  = M.r2;
  Results{i,1}.e2  = M.e2;
  Results{i,1}.rho = M.rho;
  Results{i,1}.sse = M.sse;

  Results{i,1}.r22  = M2.r2;
  Results{i,1}.e22  = M2.e2;
  Results{i,1}.rho2 = M2.rho;
  Results{i,1}.sse2 = M2.sse;  
  Results{i,1}.logModelPosteriors = F.logModelPosteriors;
end

if (isa(logModelPriors,'function_handle'))
  outputFileName = ['resultsBestModel3Bayes1Laplace' sprintf('_slag%03d_sreg%03d_maxLag%d_funcModelPrior.mat', round(sigmaLag*100), round(sigmaReg*100), lagLimit)];
else
  outputFileName = ['resultsBestModel3Bayes1Laplace' sprintf('_slag%03d_sreg%03d_maxLag%d_constPrior%d.mat', round(sigmaLag*100), round(sigmaReg*100), lagLimit, logModelPriors(1))];
end

outputFileName = fullfile(targetDir, outputFileName);
save(outputFileName,'Results', 'fitOptions', 'dataFile');

function outputFileName = ComputeColormaps(modelFitResultsFile)
global targetDir

Output      = load(modelFitResultsFile);
[r2,e2,snr] = ComputeLinearityIndicesFromBayes1LaplaceResults2(Output.Results, Output.dataFile);

r2 = reshape(r2, 168, []);
e2 = reshape(e2, 168, []);
snr= reshape(snr,168, []);

disp('Plotting the clickmap...');
[colorMaps, I] = PlotBayes1LaplaceResultsAsClickMap(Output.Results, Output.dataFile,'fitOptions',Output.fitOptions,'e2',e2,'r2',r2,'snrDb',10*log10(snr),'plotConcSeries',false);

outputFileName = fullfile(targetDir, 'dataForFigure');
save(outputFileName, 'colorMaps', 'I');
disp('Done.');

function outputFileName = ComputeSnrAngleData(modelFitResultsFile)
global targetDir;

Output = load(modelFitResultsFile);

morphInds = 1:11;
R = reshape(Output.Results, 168, []);
R = R(:,morphInds);
bm = cellfun(@(x) x.bestModel, R);
bm = bm(:);
modelTypes = {[2,9,5,12],... % octanal (unit, scaled), + lagged
              [3,10,6,13],... % citral  (unit, scaled), + lagged
              [4,11,7,14,8,15]}; % mixture  (unit, scaled, free), + lagged

ct = 0*bm;
for i = 1:numel(modelTypes)
  ct(find(sum(bsxfun(@eq, bm, modelTypes{i}),2))) = i;
end

[r2,foo,snr,signalPower, noisePower] = ComputeLinearityIndicesFromBayes1LaplaceResults2(Output.Results, Output.dataFile);

meanNoisePower = mean(noisePower,2);
snr0    = snr;
snr     = bsxfun(@rdivide, signalPower, meanNoisePower);
snrDb   = 10*log10(snr);
indResp = find(r2(:)>0.02);

th = atan2(snrDb(:,2), snrDb(:,1));
ma = sqrt(sum(snrDb(:,1:2).^2,2));
th = th(1:168*numel(morphInds));
  
snrDbResp = snrDb(:,end);
snrDbResp = snrDbResp(1:168*numel(morphInds));

outputFileName = fullfile(targetDir, 'snrAngleData.mat');
save(outputFileName, 'th', 'ct');

