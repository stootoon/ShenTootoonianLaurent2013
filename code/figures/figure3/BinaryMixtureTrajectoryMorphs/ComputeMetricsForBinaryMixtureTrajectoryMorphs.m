function Metrics = ComputeMetricsForBinaryMixtureTrajectoryMorphs(t0, binSize, responseWindow, numBs, yrDist, yrPAF, coords, nmc)
% Metrics = ComputeMetricsForBinaryMixtureTrajectoryMorphs(t0, binSize, responseWindow, numBs, yrDist, yrPAF, coords)
%
% In this function we compute the binary mixture morphs metrics. There are four metrics:
% 
% 1. Normalized Euclidean Distance (global)
% 1. Normalized Correlation Distance (global)
% 2. PAF to Citral (global)
% 3. Normalized Euclidean Distance (per-bin)
% 3. Normalized Correlation Distance (per-bin)
% 4. PAF to Citral (per-bin)

mixtureLevels = [0   140; 30  140; 60  140; 80  140; 100 140; 120 140; 140 140;
                 140 120; 140 100; 140 80;  140  60; 140  30; 140   0];
                 
mixtureInds = arrayfun(@(i) GetIndexForBinaryMixtureConcentrationPair(mixtureLevels(i,1),mixtureLevels(i,2)), 1:13);
numMixtures = numel(mixtureInds);
numNonPure = numel(mixtureInds(2:12));

pcCoords  = mixtureLevels(2:12,1)./sum(mixtureLevels(2:12,:),2);
logCoords = log10(mixtureLevels(2:12,1)./mixtureLevels(2:12,2)); 

coordinateName = when(isequal(coords,'log'),'LOG','PERCENTAGE');
x = when(isequal(coords, 'log'), logCoords, pcCoords);
xr = [min(x) max(x)];

posteriorFuncs = {@(x,y,xr,yr) ComputePosteriorOfConstantModelDirectly(x,y,xr,yr),...
                  @(x,y,xr,yr) ComputePosteriorOfLinearModel(x,y,xr,yr),...
                  @(x,y,xr,yr) ComputePosteriorOfSingleStepModelDirectly(x,y,xr,yr),...
                  @(x,y,xr,yr) ComputePosteriorOfTwoStepModel(x,y,xr,yr,nmc)};

%% Compute the distances and PAF
% Compute the Parameters structure we will use.

if (t0>responseWindow(1))
  error('Starttime is greater than the start of the response window.');
end
t1 = responseWindow(2)+2*binSize;

fprintf('Computing Parameters structure...\n');
Parameters = ComputeInputParametersForTrajectoryComputations_noduplicates('BinaryMixture','startTime',t0,'endTime',t1,'binSize',binSize,'sampleType', numBs);
fprintf('Done.\n');

Isample = ComputeTrialSamplingFromParametersStructure(Parameters);
numSamples = numel(Isample);

t = Parameters.startTime:Parameters.binSize:Parameters.endTime;
responseBins = find(t>=responseWindow(1) & t<responseWindow(2));
numResponseBins = numel(responseBins);

fprintf('Computing trajectory distances...');
[DeucOv,  DeucPb]  = ComputeTrajectoryDistancesToBinaryMixture(Parameters, [0 140], responseBins, responseBins);
[DcorrOv, DcorrPb] = ComputeTrajectoryDistancesToBinaryMixture(Parameters, [0 140], responseBins, responseBins, 'correlation');
fprintf('done.\n');

fprintf('Computing projection angle fractions...');
[PAFov, PAFpb] = ComputePafOverallAndPerBin(Parameters, responseBins, mixtureInds);
fprintf('done.\n');

M = struct;

%% Metric 1: Fit of normalized distance between response window trajectories.
% Normalize the distances.
D = bsxfun(@minus,   DeucOv, DeucOv(1,:));
D = bsxfun(@rdivide, D, D(end,:));
D = D(2:12,:);
% Compute the metric
yr = yrDist;
bestModel = zeros(1, numSamples);
logp = zeros(numSamples,4);
Y = zeros(numSamples, size(D,1));
parfor i = 1:numSamples
  y = D(:,i);
  lp = cellfun(@(f) f(x,y,xr,yr), posteriorFuncs);
  logp(i,:) = lp;
  bestModel(i) = find(lp==max(lp));
  Y(i,:) = y;
end
fprintf('The %s of [x] vs OVERALL euclidean distance was fit best by a LINEAR function in %1.0f%% of trial samplings.\n', coordinateName, mean(bestModel==2)*100);
M(1).data = 'Overall euclidean distance from citral.';
M(1).logp = logp;
M(1).bestModel = bestModel;
M(1).x = x;
M(1).Y = Y;
M(1).yr = yr;

%% Metric 2: Fit of normalized correlation distance between response window trajectories.
% Don't normalize the distances.
D = DcorrOv(2:12,:);
% Compute the metric
yr = [0 2];
bestModel = zeros(1, numSamples);
logp = zeros(numSamples,4);
Y = zeros(numSamples, size(D,1));
parfor i = 1:numSamples
  y = D(:,i);
  lp = cellfun(@(f) f(x,y,xr,yr), posteriorFuncs);
  logp(i,:) = lp;
  bestModel(i) = find(lp==max(lp));
  Y(i,:) = y;
end
fprintf('The %s of [x] vs OVERALL correlation distance was fit best by a LINEAR function in %1.0f%% of the trial samplings.\n', coordinateName, mean(bestModel==2)*100);
M(2).data = 'Overall correlation distance from citral.';
M(2).logp = logp;
M(2).bestModel = bestModel;
M(2).x = x;
M(2).Y = Y;
M(2).yr = yr;

%% Metric 3: Fit of PaF between response window trajectories.
bestModel = zeros(1, numSamples);
logp = zeros(numSamples,4);
yr = yrPAF;
Y = zeros(numSamples, size(PAFov,1));
parfor i = 1:numSamples
  y = PAFov(:,i);
  lp = cellfun(@(f) f(x,y,xr,yr), posteriorFuncs);
  logp(i,:) = lp;
  bestModel(i) = find(lp==max(lp));
  Y(i,:) = y;
end
fprintf('The %s of [x] vs OVERALL PAF was fit best by a LINEAR function in %1.0f%% of the trial samplings.\n', coordinateName, mean(bestModel==2)*100);
M(3).data = 'Overall PaF from citral.';
M(3).logp = logp;
M(3).bestModel = bestModel;
M(3).x = x;
M(3).Y = Y;
M(3).yr = yr;

%% Metric 4: Compute the pointwise fits of the euclidean in the
% response window and report the fraction of the window for which the
% LINEAR model was the best.

% Normalize the distances.
D = reshape(DeucPb,numMixtures,[]);
D = bsxfun(@minus, D, D(1,:));
D = bsxfun(@rdivide, D, D(end,:));
D = reshape(D, numMixtures, numResponseBins, numSamples);

% Compute the metric.
yr = yrDist;
bestModel = zeros(numSamples, numResponseBins);
logp      = zeros(numSamples, numResponseBins, 4);
Y = zeros(numSamples, numResponseBins, 11); 
for i = 1:numSamples
  ProgressDot(i,numSamples,25);
  YY = zeros(numResponseBins, 11);
  parfor j = 1:numResponseBins      
    y = squeeze(D(2:12, j, i));
    lp = cellfun(@(f) f(x,y,xr,yr), posteriorFuncs);
    logp(i,j,:) = lp;
    bestModel(i,j) = find(lp==max(lp));
    YY(j,:) = y;
  end
  Y(i,:,:) = YY;
end
M(4).data = 'Per bin euclidean distance from citral.';
M(4).logp = logp;
M(4).bestModel = bestModel;
M(4).x = x;
M(4).Y = Y;
M(4).yr = yr;


%% Metric 5: Compute the pointwise fits of the correlation distance in the
% response window and report the fraction of the window for which the
% LINEAR model was the best.

% Don't Normalize the distances.
D = DcorrPb;
% Compute the metric.
yr = [0 2];
bestModel = zeros(numSamples, numResponseBins);
logp      = zeros(numSamples, numResponseBins, 4);
Y = zeros(numSamples, numResponseBins, numNonPure);
for i = 1:numSamples
  ProgressDot(i,numSamples,25);
  YY = zeros(numResponseBins, numNonPure);
  parfor j = 1:numResponseBins      
    y = squeeze(D(2:12, j, i));
    lp = cellfun(@(f) f(x,y,xr,yr), posteriorFuncs);
    logp(i,j,:) = lp;
    bestModel(i,j) = find(lp==max(lp));
    YY(j,:) = y;
  end
  Y(i,:,:) = YY;
end
M(5).data = 'Per bin correlation distance from citral.';
M(5).logp = logp;
M(5).bestModel = bestModel;
M(5).x = x;
M(5).Y = Y;
M(5).yr = yr;

%% Metric 6: Compute the pointwise fits of concentration to PAF in the
% response window and report the fraction of the window for which the 
% LINEAR model was the best.

yr = yrPAF;
bestModel = zeros(numSamples, numResponseBins);
logp      = zeros(numSamples, numResponseBins, 4);
Y = zeros(numSamples, numResponseBins, numNonPure);
for i = 1:numSamples
  ProgressDot(i,numSamples,25);
  YY = zeros(numResponseBins, numNonPure);
  parfor j = 1:numResponseBins      
    y = squeeze(PAFpb(:, j, i));
    lp = cellfun(@(f) f(x,y,xr,yr), posteriorFuncs);
    logp(i,j,:) = lp;
    bestModel(i,j) = find(lp==max(lp));
    YY(j,:) = y;
  end
  Y(i,:,:) = YY;
end

M(6).data = 'Per bin PaF to citral.';
M(6).logp = logp;
M(6).bestModel = bestModel;
M(6).x = x;
M(6).Y = Y;
M(6).yr = yr;

%% Packup the results in a structure and return
Metrics.t0 = t0;
Metrics.binSize = binSize;
Metrics.yrDist = yrDist;
Metrics.yrPAF = yrPAF;
Metrics.coords = coords;
Metrics.nmc = nmc;
Metrics.x = x;
Metrics.xr = xr;
Metrics.Results = M;
Metrics.computedDate = datestr(now);
Metrics.Ibs = Parameters.Ibs;