function MakeMetricsForPaper(varargin)
% function MakeMetricsForPaper(varargin)

p = inputParser;
p.addOptional('dataDir', 'originalData');
p.parse(varargin{:});

figDir  = GetDataDirForFigure(4);
currDir = ''; % No subdirectories for this figure.
dataDir = p.Results.dataDir;

sourceDir = fullfile(figDir, currDir, dataDir);

LOADF  = @(fileName) load(fullfile(sourceDir, fileName));

MUSDSE = @(x) [mean(x) std(x) std(x)/sqrt(numel(x))];

resultsFile = 'bestModelsForMixturesResults_allOdors_sigmaReg_1_0_sigmaLag_1_0_lagLimit_3_vAlpha_1_0_vBeta_1_0.mat';
Output = LOADF(resultsFile);

mc = cellfun(@(x) sum(x{2}),   Output.Results);
r2 = cellfun(@(x) x{3}(2),     Output.Results);
ml = cellfun(@(x) numel(x{2}), Output.Results);
[snr, SignalPower, NoisePower] = ComputeSnrForAllCellsAndMixtures('dataFile',fullfile(sourceDir, 'pnComplexMixtures0_0to13_0_50ms.mat'));

%% What fraction of the models are fit?
notConst = mc~=0;
fracFit = mean(notConst(:));
fracFitPerMl = accumarray(ml(:), notConst(:), [], @mean);
fracFitPerMl = fracFitPerMl(2:end);
fracOdorAcFit = mean(notConst(:,2));
fprintf('Fraction of responses fit: %1.3f\n', fracFit);
fprintf('Fraction of responses fit at MLs 2,3,4,5 and 8: %1.3f, %1.3f, %1.3f, %1.3f, %1.3f.\n', fracFitPerMl(1), fracFitPerMl(2), fracFitPerMl(3), fracFitPerMl(4), fracFitPerMl(end));
fprintf('Fraction of odor AC fit: %1.3f\n', fracOdorAcFit);

%% Fraction fit above 3 dB
snr = snr(:,1:size(mc,2));
snrDb = 10*log10(snr);
above3dB = snrDb>3;
fracAbove3dB = mean(above3dB(:));
fracFitAbove3dB = mean(notConst(above3dB));
fracOdorAcFitAbove3dB = mean(notConst(find(above3dB(:, 2)),2));
fprintf('Fraction of responses above 3dB: %1.3f\n', fracAbove3dB);
fprintf('Fraction fit above 3dB:       %1.3f\n', fracFitAbove3dB);
fprintf('Fraction odor AC fit above 3dB: %1.3f\n', fracOdorAcFitAbove3dB);

%% Variance explained
r2 = r2(:,1:size(mc,2));
r2(isinf(r2)) = min(r2(~isinf(r2)));

muR2 = mean(r2(:));
sdR2 = std(r2(:));
seR2 = sdR2/sqrt(numel(r2));

muR2above3dB = mean(r2(above3dB));
sdR2above3dB = std(r2(above3dB));
seR2above3dB = sdR2/sqrt(numel(r2(above3dB)));

muR2notConst = mean(r2(notConst));
sdR2notConst = std(r2(notConst));
seR2notConst = sdR2/sqrt(numel(r2(notConst)));

fprintf('Mean +/- SEM of R2 overall:   (%1.2e, %1.2e)\n', muR2, seR2);
fprintf('Mean +/- SEM of R2 above 3dB: (%1.2e, %1.2e)\n', muR2above3dB, seR2above3dB);
fprintf('Mean +/- SEM of R2 for fits:  (%1.2e, %1.2e)\n', muR2notConst, seR2notConst);

%% Model distributions
ModelStats = LoadVarFromMatFileByName(fullfile(sourceDir, 'modelStats.mat'),'ModelStats');

modelTypes = zeros(8,4); % 8 components, 4 types: unit, unit-lagged, scaled, scaled-lagged
ilag = [ModelStats.ilag];
itype= [ModelStats.itype];
istruct= [ModelStats.istruct];

indSingle = find(itype>=1 & itype<9);
subs = [itype(indSingle)' istruct(indSingle)'+2*(ilag(indSingle)>0)'];
singleModelDistribution = accumarray(subs,1);
singleModelDistribution = singleModelDistribution(:,[1 3 2 4]); % 1,1-lag, k,k-lag

% Compute some metrics about the model distribution
fracUnit = sum(sum(singleModelDistribution(:,[1 2])))/sum(singleModelDistribution(:));
fracLag = sum(sum(singleModelDistribution(:,[2 4])))/sum(singleModelDistribution(:));
fprintf('Fraction of unit (and not) models: %1.3e, %1.3e\n',   fracUnit, 1-fracUnit);
fprintf('Fraction of lagged (and not) models: %1.3e, %1.3e\n', fracLag,  1-fracLag ); 

%% Scaling weights
ModelStats = LoadVarFromMatFileByName(fullfile(sourceDir, 'modelStats.mat'),'ModelStats');
dataForFigure = LOADF('dataForFigure.mat');

B = GetOdorNamesAsBinaryVectors('full');
B = B(13:end,:);
ml = sum(B,2);
  
istruct = [ModelStats.istruct];
istruct = reshape(istruct, size(ModelStats));
istructIsScaled = (istruct == 1 | istruct==2); % Look for the scaled models
counts = zeros(8,4); % 8 components, 4 differences (3-2,4-2,5-2,8-2)
diffs  = zeros(8,4);
mls = [2:5 8];
recs = [];

for i = 1:8 % Loop over components
  inds = dataForFigure.I1{i};
  for j = 1:numel(inds) % Loop over cells that produced responses of this type
    newRec = nan(1,5);
    for iml = 1:5    % Find the mean k at each ml
      thisMl = mls(iml);
      indMl = find(ml==thisMl & istructIsScaled(inds(j),:)'); % Get all the responses that were scaled or unit.
      if (~isempty(indMl)) % Found some at this ml
        w1 = zeros(1,numel(indMl));
        for k = 1:numel(indMl) % Grab the weights
          w = ModelStats(inds(j),indMl(k)).weights;
          w = unique(w(~isnan(w)));
          assert(~isempty(w) && numel(w)==1,'Expected exactly 1 weight for scaled model.');
          w1(k) = w;
        end
        newRec(iml) = mean(w1); % Store their mean value
      end
    end
    indNonNan = find(~isnan(newRec));    
    recs(size(recs,1)+1,:) = [i j newRec];
  end
end
fprintf('%d PNs used.\n', size(recs,1));
muK = nanmean(Columnize(recs(:,3:end)));
sdK = nanstd(Columnize(recs(:,3:end)));
seK = sdK/sqrt(sum(~isnan(recs(:))));
fprintf('Mean +/- SEM of scaling factors: %1.3e, %1.3e.\n', muK, seK);

% Now prune to those that have at least three ml responses
vldCount = sum(~isnan(recs(:,end-4:end)),2);
indToUse = find(vldCount>2);
recsToUse = recs(indToUse,:);
corrCoef = 0*recsToUse(:,1);
fprintf('%d/%d records used.\n', size(recsToUse,1), size(recs,1));

for i = 1:size(recsToUse,1)
  indNonNan = find(~isnan(recsToUse(i,end-4:end)));
  x = mls(indNonNan)';
  y = recsToUse(i,indNonNan+2)';
  corrCoef(i) = corr(x, y,'type','spearman');
end
fprintf('%d/%d corr coefs are nan.\n', sum(isnan(corrCoef)), numel(corrCoef));
corrCoef(isnan(corrCoef)) = 0;
fprintf('Mean Spearman Rho: %1.3e\n',   mean(corrCoef));
fprintf('Median Spearman Rho: %1.3e\n', median(corrCoef));
fprintf('P-value (signtest): %1.3e\n', signtest(corrCoef));

%% Metrics about snr angles
data = LOADF('snrAngleData.mat');
fprintf('Overall fraction of angles within pi/8 of zero (of %d total): %1.3e\n', numel(data.thAll), mean(abs(data.thAll)<pi/8));
fprintf('Fraction of angles within pi/8 of zero (preferred component responses only, %d total): %1.3e\n', numel(data.thPref), mean(abs(data.thPref)<pi/8));

%% Correlations based on weight ratio
data = LOADF('weightRatioData.mat');

corrCoef    = data.corrCoef;
corrCoef2   = data.corrCoef2;
weightRatio = data.weightRatio;
numSingleCmpFits = sum(data.isingle(:));

cutoff = 0.2;

corrBelow = corrCoef(weightRatio<cutoff);
corrAbove = corrCoef(weightRatio>=cutoff);

corr2Below = corrCoef2(weightRatio<cutoff);
corr2Above = corrCoef2(weightRatio>=cutoff);

musdseBelow = MUSDSE(corrBelow);
musdseAbove = MUSDSE(corrAbove);

musdseBelow2 = MUSDSE(corr2Below);
musdseAbove2 = MUSDSE(corr2Above);

fprintf('Fraction with corrcoef above cutoff: %1.3e (%1.3e of all single component fits)\n', numel(corrAbove)/numel(corrCoef), numel(corrAbove)/numSingleCmpFits);
fprintf('Fraction with corrcoef below cutoff: %1.3e (%1.3e of all single component fits)\n', numel(corrBelow)/numel(corrCoef), numel(corrBelow)/numSingleCmpFits);

fprintf('Mean +/- SEM of corrcoef below: %1.3e + %1.3e \n', musdseBelow(1), musdseBelow(3));
fprintf('Mean +/- SEM of corrcoef above: %1.3e + %1.3e \n', musdseAbove(1), musdseAbove(3));

fprintf('Mean +/- SEM of corr2coef below: %1.3e + %1.3e \n', musdseBelow2(1), musdseBelow2(3));
fprintf('Mean +/- SEM of corr2coef above: %1.3e + %1.3e \n', musdseAbove2(1), musdseAbove2(3));

end
