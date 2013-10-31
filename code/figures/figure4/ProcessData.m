function ProcessData()
% function ProcessData()

global targetDir

figDir    = GetDataDirForFigure(4);
currDir   = ''; % No sub directories for this figure.
targetDir = fullfile(figDir, currDir, 'recomputedData');

t0      = 0;
t1      = 6;
t1long  = 13;
binSize = 0.050;

numMpCores = 6;

fprintf('Preparing Cbot data file in (%d-%d) interval...\n', t0, t1); tic;
pnCbotDataFile = PreparePnCbotDatasets(t0, t1,     binSize);
fprintf('Done in %1.1f seconds.\n', toc);

fprintf('Computing model fits...\n'); tic;
modelFitResultsFile = ParComputeBestModelLaplaceGeneral(numMpCores, pnCbotDataFile); % ~71000 secs
fprintf('Done in %1.1f seconds.\n', toc);

fprintf('Computing model stats...\n'); tic;
modelStatsFile     = ComputeModelStats(pnCbotDataFile, modelFitResultsFile); % ~130 seconds
fprintf('Done in %1.1f seconds.\n', toc);

fprintf('Preparing Cbot data file in (%d-%d) interval...\n', t0, t1long); tic;
pnCbotDataFileLong = PreparePnCbotDatasets(t0, t1long, binSize); % ~30 secs. Will need this for the SNR computations
fprintf('Done in %1.1f seconds.\n', toc);

fprintf('Computing data for figure...\n'); tic;
dataForFigureFile  = ComputeDataForFigure(modelFitResultsFile, pnCbotDataFileLong, pnCbotDataFile); % ~4.5 seconds
fprintf('Done in %1.1f seconds.\n', toc);

fprintf('Computing SNR angle data.\n'); tic;
snrAngleDataFile   = ComputeSnrAngleData(pnCbotDataFileLong, modelStatsFile, dataForFigureFile); % ~19.5 seconds
fprintf('Done in %1.1f seconds.\n', toc);

fprintf('Preparing Weight Ratio data.\n'); tic;
weightRatioDataFile= ComputeWeightRatioData(pnCbotDataFile, modelFitResultsFile, modelStatsFile, dataForFigureFile, pnCbotDataFileLong); % ~130 secs
fprintf('Done in %1.1f seconds.\n', toc);

function outputFile = PreparePnCbotDatasets(t0, t1, binSize)
global targetDir

pnSpt   = LoadTocSpikeTimes('rawpn');
t0fit   = 1.5;
t1fit   = 5;

numOdors = 44;
numPns   = 174;

pnCbot = CountSpikesInBinsAndAverageAcrossTrials(pnSpt, arrayfun(@Identity,1:7,'UniformOutput', false), 1:numOdors, 1:numPns, 'binSize', binSize, 'startTime', t0, 'endTime', t1);

[numCells,numBins,numTrials,numOdors] = size(pnCbot);

binStarts = (0:numBins-1)*binSize;

iall = 1:numBins;
ifit = find(binStarts>=t0fit & binStarts<t1fit);
tfit = binStarts(ifit);
tall = binStarts;

datasetName = sprintf('pnComplexMixtures%1.1fto%1.1f_%dms', t0, t1, round(binSize*1000));
datasetName = strrep(datasetName, '.','_');

outputFile = fullfile(targetDir, [datasetName '.mat']);
save(outputFile, 'binSize','t0','t1','t0fit','t1fit','iall','ifit','tfit','tall', 'pnCbot');

function outputFile = ParComputeBestModelLaplaceGeneral(numMpCores, dataFile)
global targetDir

% dataFile = 'pnComplexMixtures0_0to6_0_50ms.mat';
Data = load(dataFile);

whichOdorIds = [13:44];
  
whichCells = [1:174];
[odorComp, odorNames] = GetOdorNamesAsBinaryVectors('full');

Results = {};

sigmaReg           = 1;
sigmaLag           = 1;
lagLimit           = 3;
variancePriorAlpha = 1;
variancePriorBeta  = 1;

binOffset      = find(Data.iall==Data.ifit(1),1);
whichBinsToFit = binOffset+(0:length(Data.ifit)-1);

% It's important to order the parameters in this way rather than the
% cells x odors because high component odors take longer to fit, and
% we want those fits to be distributed evenly between the cores.
params = CartesianProduct({1:numel(whichOdorIds), 1:numel(whichCells)});

numParams = size(params,1);
% Start the computation
if (MatlabPoolWrapper('size')>0)
  MatlabPoolWrapper close force local
end
MatlabPoolWrapper('open',numMpCores);

logFiles = {};
for i = 1:numMpCores
  logFile{i} = sprintf('worker%d.log',i);
  fp = fopen(logFile{i},'w');
  fclose(fp);
end

Results = cell(numParams, 1);
parfor i = 1:numParams
  currTask = getCurrentTask();
  workerId = get(currTask, 'ID');
    
  odorId    = whichOdorIds(params(i,1));  
  whichCell = whichCells(  params(i,2));

  fp = fopen(logFile{workerId}, 'a');
  fprintf(fp, 'Started fit of cell %d to odor %d.\n', whichCell, odorId);
  fclose(fp);
  
  indCmps = find(odorComp(odorId,:))+3; % Skip the high odors
  numCmps = numel(indCmps);
  
  Xobs = squeeze(Data.pnCbot(whichCell,Data.iall,indCmps,:));
  yobs = squeeze(Data.pnCbot(whichCell,Data.iall,odorId,:));
  
  Fit = FitAllModelsForObservations2LaplaceGeneral(Xobs, yobs, 'whichBinsToFit', whichBinsToFit, 'sigmaReg', sigmaReg, 'sigmaLag', sigmaLag, 'lagLimit', lagLimit, 'variancePriorAlpha',variancePriorAlpha, 'variancePriorBeta', variancePriorBeta);

  fitSummary = struct;
  fitSummary.inputConfigs   = Fit.inputConfigs;
  fitSummary.logLikelihoods = Fit.logLikelihoods;
  fitSummary.logPosteriors  = Fit.logPosteriors;
  fitSummary.exitFlags      = Fit.exitFlags;
  
  if (any(Fit.exitFlags ~= 1 & ~isinf(Fit.exitFlags)))
    fp = fopen(logFile{workerId},'a');
    fprintf(fp, 'Fit of cell %d to odor %d produced some incomplete fits.\nExit flags: ', whichCell, odorId);
    for j1 = 1:size(Fit.exitFlags,1)
      for j2 = 1:size(Fit.exitFlags,2)
        fprintf(fp,'%d ', Fit.exitFlags(j1,j2));
      end
      fprintf(fp,'\n');
    end
    fprintf(fp,'\n\n');
    fclose(fp);
  end
  
  [imax,jmax] = find(Fit.logPosteriors == max(Fit.logPosteriors(:)), 1);
  bestInputConfig = Fit.inputConfigs(imax,1:end-1);
  bestFitModel    = jmax;
  bestFitResults  = [Fit.bestModel.sse, Fit.bestModel.r2, Fit.bestModel.e2];
  Results{i} = {bestFitModel, bestInputConfig, bestFitResults, odorId, whichCell, fitSummary};
end
Results = reshape(Results, numel(whichOdorIds), numel(whichCells));
Results = Results';

% Save the outputs
odorAllStr     = when(max(whichOdorIds)==44, 'allOdors', 'noOdorALL');
saveFileName   = sprintf('bestModelsForMixturesResults_%s_sigmaReg_%1.1f_sigmaLag_%1.1f_lagLimit_%d_vAlpha_%1.1f_vBeta_%1.1f', odorAllStr, sigmaReg, sigmaLag, lagLimit, variancePriorAlpha,variancePriorBeta);
saveFileName   = strrep(saveFileName,'.','_');

outputFile = fullfile(targetDir, [saveFileName '.mat']);
save(outputFile, 'Results', 'whichOdorIds', 'whichCells', 'sigmaReg', 'sigmaLag', 'lagLimit', 'whichBinsToFit', 'variancePriorAlpha', 'variancePriorBeta', 'dataFile');

function outputFile = ComputeModelStats(pnCbotDataFile, modelFitResultsFile)
% In this script we use the bestModels results to go and actually
% compute all the weights and other stats for each model.

% pnCbotDataFile = 'pnComplexMixtures0_0to6_0_50ms.mat';
global targetDir

Data = load(pnCbotDataFile);
whichBinsToFit = find(Data.iall==Data.ifit(1))+(0:numel(Data.ifit)-1);

% resultsFile = 'bestModelsForMixturesResults_allOdors_sigmaReg_1_0_sigmaLag_1_0_lagLimit_3_vAlpha_1_0_vBeta_1_0_noduplicates.mat';
Rs = load(modelFitResultsFile);
fitOptions = {'sigmaReg',Rs.sigmaReg, 'sigmaLag', Rs.sigmaLag, 'lagLimit', Rs.lagLimit, 'variancePriorAlpha', Rs.variancePriorAlpha, 'variancePriorBeta', Rs.variancePriorBeta};

ModelStats = [];

% We have to watch out in classifying models, because 1 component
% scaled models are equivalent to free models. So we'll take all 1
% component free models to be scaled, and let the rest be whatever
% they are.
%
cmps = 'ABCDWXYZ';
for i = 1:size(Rs.Results,1)
  ProgressDot(i,size(Rs.Results,1),25);
  for j = 1:size(Rs.Results,2)
    thisOdorCmps = GetComponentsForOdor(j+12);
    whichInputConfig = Rs.Results{i,j}{2};
    numInputsUsed    = sum(whichInputConfig);
    whichInputsUsed  = find(whichInputConfig);
    
    if (numInputsUsed > 0 ) % A non constant model 
      Xobs = squeeze(Data.pnCbot(i, :, thisOdorCmps+3, :));
      yobs = squeeze(Data.pnCbot(i, :, j+12, :));

      whichModel = Rs.Results{i,j}{1};
      FitBest = FitAllModelsForObservations2LaplaceGeneral(Xobs, yobs, 'whichBinsToFit', whichBinsToFit, fitOptions{:}, 'fitMode', 'single', 'whichInputConfig',whichInputConfig,'whichModel',whichModel);           
      
      w = nan(1,8);      
      if (numInputsUsed == 1) % A single-component model
        switch(whichModel)
         case {1,4} % Unit Model
          w(thisOdorCmps(whichInputsUsed)) = 1;
          ModelStats(i,j).structure = 'unit';
          ModelStats(i,j).istruct = 1;
          ModelStats(i,j).type = cmps(thisOdorCmps(whichInputsUsed));
          ModelStats(i,j).itype= thisOdorCmps(whichInputsUsed);
          ModelStats(i,j).lagged = when(whichModel<=3,'no','yes');
          ModelStats(i,j).ilag   = double(whichModel>3);
          ModelStats(i,j).weights = w;
         case {2,3,5,6} % Scaled Model (Combined with free)
          w(thisOdorCmps(whichInputsUsed)) = FitBest.Model.w(end);
          ModelStats(i,j).structure = 'scaled';
          ModelStats(i,j).istruct = 2;
          ModelStats(i,j).type = cmps(thisOdorCmps(whichInputsUsed));
          ModelStats(i,j).itype= thisOdorCmps(whichInputsUsed);
          ModelStats(i,j).lagged = when(whichModel<=3,'no','yes');
          ModelStats(i,j).ilag   = double(whichModel>3);
          ModelStats(i,j).weights = w;
         otherwise
          error('Unknowned model type: %d\n', whichModel);
        end
        ModelStats(i,j).numInputsUsed = 1;
      else % A multi-component model
        switch(whichModel)
         case {1,4}
          w(thisOdorCmps(whichInputsUsed)) = 1;
          ModelStats(i,j).structure = 'unit';
          ModelStats(i,j).istruct = 1;
          ModelStats(i,j).type = cmps(thisOdorCmps(whichInputsUsed));
          ModelStats(i,j).itype= 10;
          ModelStats(i,j).lagged = when(whichModel<=3,'no','yes');
          ModelStats(i,j).ilag   = double(whichModel>3);
          ModelStats(i,j).weights = w;          
         case {2,5}
          w(thisOdorCmps(whichInputsUsed)) = FitBest.Model.w(end);
          ModelStats(i,j).structure = 'scaled';
          ModelStats(i,j).istruct = 2;
          ModelStats(i,j).type = cmps(thisOdorCmps(whichInputsUsed));
          ModelStats(i,j).itype= 10;
          ModelStats(i,j).lagged = when(whichModel<=3,'no','yes');
          ModelStats(i,j).ilag   = double(whichModel>3);
          ModelStats(i,j).weights = w;
         case {3,6}
          w(thisOdorCmps(whichInputsUsed)) = FitBest.Model.w(2:end);
          ModelStats(i,j).structure = 'free';
          ModelStats(i,j).istruct = 3;
          ModelStats(i,j).type = cmps(thisOdorCmps(whichInputsUsed));
          ModelStats(i,j).itype= 10;
          ModelStats(i,j).lagged = when(whichModel<=3,'no','yes');
          ModelStats(i,j).ilag   = double(whichModel>3);
          ModelStats(i,j).weights = w;
        end
      end
      ModelStats(i,j).bestModel    = whichModel;
      ModelStats(i,j).inputConfigs = whichInputConfig;
      ModelStats(i,j).descr = FitBest.descr;
      ModelStats(i,j).r2 = FitBest.Model.r2;
      ModelStats(i,j).e2 = FitBest.Model.e2;
      ModelStats(i,j).rho= FitBest.Model.rho;
      ModelStats(i,j).sse= FitBest.Model.sse;
      ModelStats(i,j).var= FitBest.Model.v;      
    else
      w(thisOdorCmps) = 0;
      ModelStats(i,j).numInputsUsed = 0;
      ModelStats(i,j).weights = w;
      ModelStats(i,j).structure = 'constant';
      ModelStats(i,j).istruct = 0;
      ModelStats(i,j).type = '-';
      ModelStats(i,j).itype= 0;
      ModelStats(i,j).lagged = 'no';
      ModelStats(i,j).ilag = 0;
    end
  end
end

outputFile = fullfile(targetDir, 'modelStats.mat');
save(outputFile, 'ModelStats');

function outputFile = ComputeDataForFigure(modelFitResultsFile, pnCbotDataFileLong, pnCbotDataFile)
global targetDir

[colorMaps, colorMapNames, pnPlotOrder, I1, Imix, Iamb, Inull] = PlotResultsAsClickmap(modelFitResultsFile, pnCbotDataFileLong, pnCbotDataFile);
outputFile = fullfile(targetDir, 'dataForFigure.mat');
save(outputFile, 'colorMaps', 'colorMapNames', 'pnPlotOrder', 'I1', 'Imix', 'Iamb', 'Inull');

function outputFile = ComputeSnrAngleData(pnDataFileLong, modelStatsFile, dataForFigureFile)
% In this script we compute the SNR angle for all the single component responses,
% and then again for the responses which were to the cell's perferred response.
global targetDir

ModelStats = LoadVarFromMatFileByName(modelStatsFile,'ModelStats');
dataForFigure = load(dataForFigureFile);

[snr, SignalPower, NoisePower] = ComputeSnrForAllCellsAndMixtures('dataFile',pnDataFileLong);
muNoise = cellfun(@mean, NoisePower);

itype = reshape([ModelStats.itype], size(ModelStats));
isingle = (itype>=1 & itype<=8);

% What we want to do is loop over the single responses, and compute
% the snr of the single component, and the maximum of the other present
[icell, jmixture] = find(isingle);  
theta = nan*isingle;
for ii = 1:numel(icell)
  ProgressDot(ii,numel(icell), 25);
  i = icell(ii);
  j = jmixture(ii);
  np = NoisePower{i,j};
  sp = SignalPower{i,j};
  muNoise = mean(np);
  mixCmps = GetComponentsForOdor(j+12);
  prefCmp = ModelStats(i,j).itype;
  indPref = find(mixCmps==prefCmp);
  assert(numel(indPref)==1,'Expected the preferred component to show up exactly once.');
  indOther = setdiff(1:numel(mixCmps), indPref);
  maxOtherSignals = max(sp(indOther));
  prefSnrDb = 10*log10(sp(indPref)/muNoise);
  maxOtherSnrDb = 10*log10(maxOtherSignals/muNoise);
  snrAngle = atan2(maxOtherSnrDb, prefSnrDb);
  theta(i,j) = snrAngle;
end
thAll = theta(~isnan(theta));

% Now loop over the cells that have a preferred response, and look
% only at their preferred responses.

theta2 = nan*theta;
for prefCmp = 1:8
  for ii = 1:numel(dataForFigure.I1{prefCmp})
    i = dataForFigure.I1{prefCmp}(ii);
    indPrefResp = find(itype(i,:)==prefCmp);
    for jj = 1:numel(indPrefResp)
      j = indPrefResp(jj);
      np = NoisePower{i,j};
      sp = SignalPower{i,j};
      muNoise = mean(np);
      mixCmps = GetComponentsForOdor(j+12);
      indPref = find(mixCmps==prefCmp);
      assert(numel(indPref)==1,'Expected the preferred component to show up exactly once.');
      indOther = setdiff(1:numel(mixCmps), indPref);
      maxOtherSignals = max(sp(indOther));
      prefSnrDb = 10*log10(sp(indPref)/muNoise);
      maxOtherSnrDb = 10*log10(maxOtherSignals/muNoise);
      snrAngle = atan2(maxOtherSnrDb, prefSnrDb);
      theta2(i,j) = snrAngle;
    end
  end
end
thPref = theta2(~isnan(theta2));

outputFile = fullfile(targetDir, 'snrAngleData.mat');
save(outputFile, 'thAll', 'thPref');

function outputFile = ComputeWeightRatioData(pnCbotDataFile, modelFitResultsFile, modelStatsFile, dataForFigureFile, pnCbotDataFileLong)
% In this script we compute the data required for Figure S2G, the
% weight ratio panel.
global targetDir

ModelStats    = LoadVarFromMatFileByName(modelStatsFile, 'ModelStats');
dataForFigure = load(dataForFigureFile);

w1 = [];
w2 = [];

[snr, SignalPower, NoisePower] = ComputeSnrForAllCellsAndMixtures('dataFile',pnCbotDataFileLong); % ComputeSnrForAllCellsAndMixtures();
muNoise = cellfun(@mean, NoisePower);

itype = reshape([ModelStats.itype], size(ModelStats));
isingle = (itype>=1 & itype<=8);

theta2 = nan*itype;
indSecondary = nan*itype; % This will store the index of the secondary model
for prefCmp = 1:8
  for ii = 1:numel(dataForFigure.I1{prefCmp})
    i = dataForFigure.I1{prefCmp}(ii);
    indPrefResp = find(itype(i,:)==prefCmp);
    for jj = 1:numel(indPrefResp)
      j = indPrefResp(jj);
      np = NoisePower{i,j};
      sp = SignalPower{i,j};
      muNoise = mean(np);
      mixCmps = GetComponentsForOdor(j+12);
      indPref = find(mixCmps==prefCmp);
      assert(numel(indPref)==1,'Expected the preferred component to show up exactly once.');
      indOther = setdiff(1:numel(mixCmps), indPref);
      maxOtherSignals = max(sp(indOther));
      indSecondary(i,j) = indOther(argmax(sp(indOther),1));
      prefSnrDb = 10*log10(sp(indPref)/muNoise);
      maxOtherSnrDb = 10*log10(maxOtherSignals/muNoise);
      snrAngle = atan2(maxOtherSnrDb, prefSnrDb);
      theta2(i,j) = snrAngle;
    end
  end
end

[irefit, jrefit] = find(abs(theta2-pi/4)<pi/8); % This used to be find(abs(theta2-pi/4)<pi/8)
fprintf('%d cases with |SNR angle - pi/4|<pi/8.\n', numel(irefit));

Data = load(pnCbotDataFile);
Rs   = load(modelFitResultsFile);

whichOdorIds = [13:44];

[odorComp, odorNames] = GetOdorNamesAsBinaryVectors('full');

corrCoef = 0*irefit;
for i = 1:numel(irefit)
  ProgressDot(i, numel(irefit), 25);
  whichCell = irefit(i);
  odorId = jrefit(i)+12;
  indCmps = find(odorComp(odorId,:))+3; % Skip the high odors
  numCmps = numel(indCmps);
  
  Xobs = squeeze(Data.pnCbot(whichCell,Data.iall,indCmps,:));
  yobs = squeeze(Data.pnCbot(whichCell,Data.iall,odorId,:));
  
  primaryConfig = ModelStats(irefit(i), jrefit(i)).inputConfigs;
  secondaryConfig = 0*primaryConfig;
  secondaryConfig(indSecondary(irefit(i), jrefit(i))) = 1;
  primaryWeight = ModelStats(irefit(i),jrefit(i)).weights;
  
  primaryWeight = primaryWeight(~isnan(primaryWeight));
  
  % Refit the primary Model
  primaryBestModel = ModelStats(irefit(i), jrefit(i)).bestModel;
  primaryFit = FitAllModelsForObservations2LaplaceGeneral(Xobs, yobs,...
                                                    'whichBinsToFit', Rs.whichBinsToFit, 'sigmaReg', Rs.sigmaReg, 'sigmaLag', Rs.sigmaLag, 'lagLimit', Rs.lagLimit, 'variancePriorAlpha',Rs.variancePriorAlpha, 'variancePriorBeta', Rs.variancePriorBeta,...
                                                    'fitMode', 'single', 'whichInputConfig',primaryConfig,'whichModel', primaryBestModel);
  if (any(primaryBestModel==[1 4]))
    w1(i) = 1;
  else
    w1(i) = primaryFit.Model.w(end);
  end
  
  maxLp = -inf;
  ibestSecondary = nan;
  Fit = {};
  for whichModel = 1:6
    Fit{whichModel} = FitAllModelsForObservations2LaplaceGeneral(Xobs, yobs,...
                                                      'whichBinsToFit', Rs.whichBinsToFit, 'sigmaReg', Rs.sigmaReg, 'sigmaLag', Rs.sigmaLag, 'lagLimit', Rs.lagLimit, 'variancePriorAlpha',Rs.variancePriorAlpha, 'variancePriorBeta', Rs.variancePriorBeta,...
                                                      'fitMode', 'single', 'whichInputConfig',secondaryConfig,'whichModel',whichModel);
    lp = Fit{whichModel}.logPosteriors;
    if (lp>maxLp)
      maxLp = lp;
      w2(i) = when(any(whichModel ==[1 4]), 1, Fit{whichModel}.Model.w(end));
      ibestSecondary = whichModel;
    end
  end  

  X1all = squeeze(mean(Xobs(:, find(primaryConfig),:),3));
  X2all = squeeze(mean(Xobs(:, find(secondaryConfig),:),3));
  
  X1fit = X1all(Rs.whichBinsToFit,:);
  X2fit = X2all(Rs.whichBinsToFit,:);
  yfit  = mean(yobs(Rs.whichBinsToFit,:),2);
  
  corrCoef(i) = corr(X1fit, X2fit);
  corrCoef2(i) = corr(X2fit, yfit);
end
weightRatio = w2./w1;

outputFile = fullfile(targetDir, 'weightRatioData.mat');
save(outputFile, 'weightRatio', 'corrCoef', 'corrCoef2', 'isingle');
