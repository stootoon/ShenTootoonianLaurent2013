function MakeMetricsForPaper(whichMetrics, varargin)
% function MakeMetricsForPaper(whichMetrics, [dataDir = originalData])
%
% whichMetrics = 1 => How many null responses?
% whichMetrics = 2 => Mean+/-SEM of R2 for non-null cells.
% whichMetrics = 3 => The distribution of scaling weights
% whichMetrics = 4 => Is there an effect of morph on the scaling weight
% whichMetrics = 5 => Is there an effect of morph on the scaling weight Part II - compare to 1, the value we'd have (by definition) at the edges 
% whichMetrics = 6 => Compute fraction fit stats at 100:100, because we can compare it to what we find for the equivalent in the complex mixtures setup, odor AC.
% whichMetrics = 7 => Fraction of responses to cit100:oct100 that are unit, and lagged
% whichMetrics = 8 => Mean +/sem of scaling weights in preferred cells for preferred responses
% whichMetrics = 10=> What fraction of models for the morph data  were fit with non-constant fits
% whichMetrics = 11=> Fraction of high strong single component responses

p = inputParser;
p.addOptional('dataDir', 'originalData');
p.parse(varargin{:});

figDir  = GetDataDirForFigure(2);
currDir = GetCurrentDirFromPathString(fileparts(mfilename('fullpath')));
dataDir = p.Results.dataDir;

sourceDir = fullfile(figDir, currDir, dataDir);

pnCbotDataFile = fullfile(sourceDir, 'delayRegressData.mat');

LOADF = @(fileName) load(fullfile(figDir, currDir, dataDir, fileName));

if (nargin==0)
  whichMetrics = 1:10;
end

%% Metric 1: How many null responses?
if (any(whichMetrics==1))
  dataForFigure = LOADF('dataForFigure.mat');
  groupSizes = cellfun(@numel, dataForFigure.I);
  fprintf('%d/%d = %1.3f %% of PNs were assigned to NULL.\n', groupSizes(1), sum(groupSizes), round(groupSizes(1)/sum(groupSizes)*100));
end

%% Metric 2: Mean+/-SEM of R2 for non-null cells.
if (any(whichMetrics==2))
  resultsFile = 'resultsBestModel3Bayes1Laplace_slag100_sreg100_maxLag3_funcModelPrior.mat'; %resultsBestModel3Bayes1Laplace_slag100_sreg100_maxLag5_constPrior3.mat';
  Output = LOADF(resultsFile);
  [r2,e2,snr,signalPower, noisePower] = ComputeLinearityIndicesFromBayes1LaplaceResults2(Output.Results, pnCbotDataFile);
  r2 = reshape(r2, 168, []);
  e2 = reshape(e2, 168, []);
  snr= reshape(snr,168, []);
  morphInds = 1:11;
  
  r2 = r2(:,morphInds);
  snr = snr(:,morphInds);
  r2 = r2(:);
  snrDb = 10*log10(snr(:));

  for snrDbTh = [-inf 0 1 2 3]
    fprintf('Mean +/- SEM of R2 when SNR>%1.2e dB: %1.2e +/- %1.2e (sd = %1.2e) \n', snrDbTh, mean(r2(snrDb>snrDbTh)), std(r2(snrDb>snrDbTh))/sqrt(sum(snrDb>snrDbTh)), std(r2(snrDb>snrDbTh)));  
  end
end

%% Metric 3: The distribution of scaling weights
if (any(whichMetrics==3))
  % Go through the results and find the distribution of scaling values.
  resultsFile = 'resultsBestModel3Bayes1Laplace_slag100_sreg100_maxLag3_funcModelPrior.mat'; %resultsBestModel3Bayes1Laplace_slag100_sreg100_maxLag5_constPrior3.mat';
  Output = LOADF(resultsFile);

  kcit = [];
  koct = [];
  kmix = [];
  morphInds = 1:11;

  R = reshape(Output.Results,168,[]);
  R = R(:,morphInds);
  for i = 1:size(R,1)
    for j = 1:size(R,2)
      switch(R{i,j}.bestModel)
       case {5,12} % Scaled Octanol
        koct(numel(koct)+1) = R{i,j}.w(2);
       case {6,13} % Scaled Citral
        kcit(numel(kcit)+1) = R{i,j}.w(2);
       case {7,14} % Scaled mixture
        kmix(numel(kmix)+1) = R{i,j}.w(2);
      end
    end
  end
  
  fprintf('Mean +/- SD, SEM scaling factor for scaled-octanol responses: %1.3f, %1.3f, %1.3e\n', mean(koct), std(koct), std(koct)/sqrt(numel(koct)));
  fprintf('Mean +/- SD, SEM scaling factor for scaled-citral  responses: %1.3f, %1.3f, %1.3e\n', mean(kcit), std(kcit), std(kcit)/sqrt(numel(kcit)));
  fprintf('Mean +/- SD, SEM  scaling factor for scaled-mixture responses: %1.3f, %1.3f, %1.3e\n', mean(kmix), std(kmix),std(kmix)/sqrt(numel(kmix)));

  fprintf('Mean +/- SD, SEM over cit and oct: %1.3f, %1.3f, %1.3e\n', mean([koct(:);kcit(:)]), std([koct(:);kcit(:)]), std([koct(:);kcit(:)])/sqrt(numel(koct)+numel(kcit)));
  fprintf('Fraction below zero or above 1, for scaled-octanol responses: %1.3f, %1.3f\n', mean(koct<0), mean(koct>1));
  fprintf('Fraction below zero or above 1, for scaled-citral  responses: %1.3f, %1.3f\n', mean(kcit<0), mean(kcit>1))
  fprintf('Fraction below zero or above 1, for scaled-mixture responses: %1.3f, %1.3f\n', mean(kmix<0), mean(kmix>1));  
end

%% Metric 4: Is there an effect of morph on the scaling weight
if (any(whichMetrics==4))
  % Go through the results and grab, for the citral-type responses
  % the k values of the scaling, and plot them.
  resultsFile = 'resultsBestModel3Bayes1Laplace_slag100_sreg100_maxLag3_funcModelPrior.mat'; %resultsBestModel3Bayes1Laplace_slag100_sreg100_maxLag5_constPrior3.mat';
  Output = LOADF(resultsFile);
  
  morphInds = 1:11;

  R = reshape(Output.Results, 168, []);
  R = R(:,morphInds);

  dataForFigure = LOADF('dataForFigure');
  
  Icit = dataForFigure.I{2};
  Ioct = dataForFigure.I{3};
  
  koct = nan(numel(Ioct), numel(morphInds));
  kcit = nan(numel(Icit), numel(morphInds));;

  R = reshape(Output.Results,168,[]);
  R = R(:,morphInds);
  
  for i = 1:numel(Ioct)
    for j = 1:size(R,2)
      switch(R{Ioct(i),j}.bestModel)
       case {2, 9} % Unit Octanol 
        koct(i,j) = 1;
       case {5,12} % Scaled Octanol
        koct(i,j) = R{Ioct(i),j}.w(2);
      end
    end  
  end
  
  muOct = zeros(1,numel(morphInds)-1);
  sdOct = zeros(1,numel(morphInds)-1);
  seOct = zeros(1,numel(morphInds)-1);
  pvOct = zeros(1,numel(morphInds)-1);
  szOct = zeros(1,numel(morphInds)-1);
  diffsOct = cell(numel(morphInds)-1,1);
  for i = 1:numel(morphInds)-1
    ind = find(~isnan(koct(:,end)) & ~isnan(koct(:,i)));
    szOct(i) = numel(ind);
    x1 = koct(ind,end);
    xn = koct(ind,i);
    diffOct{i} = xn - x1;
    muOct(i) = mean(xn-x1);
    sdOct(i) = std(xn-x1);
    seOct(i) = sdOct(i)/sqrt(szOct(i));
    [h,pvOct(i)] = ttest(x1, xn);
  end
  fprintf('Per/mixture level p-values for octanol: ');
  fprintf('%1.3f ', pvOct);
  fprintf('\n');

  % Do ANOVA to get effect of differences
  mug = sum(cellfun(@sum, diffOct))/sum(cellfun(@numel, diffOct));
  sst = sum(cellfun(@(x) sum((x-mug).^2), diffOct));
  sse = sum(cellfun(@(x) sum((x - mean(x)).^2), diffOct));
  ssr = sum(cellfun(@(x) numel(x)*(mean(x)-mug).^2, diffOct));
  dfe = sum(cellfun(@numel, diffOct)) - numel(diffOct);
  dfr = numel(diffOct) - 1;
  msr = ssr/dfr;
  mse = sse/dfe;
  fstat = msr/mse;
  pvOctOverall = 1 - fcdf(fstat, dfr, dfe);
  
  for i = 1:numel(Icit)
    for j = 1:size(R,2)
      switch(R{Icit(i),j}.bestModel)
       case {3, 10} % Unit Citral 
        kcit(i,j) = 1;
       case {6, 13} % Scaled Citral
        kcit(i,j) = R{Icit(i),j}.w(2);
      end
    end
  end

  muCit = zeros(1,numel(morphInds)-1);
  sdCit = zeros(1,numel(morphInds)-1);
  seCit = zeros(1,numel(morphInds)-1);
  pvCit = zeros(1,numel(morphInds)-1);
  szCit = zeros(1,numel(morphInds)-1);
  diffsCit = cell(numel(morphInds)-1,1);   
  for i = 1:numel(morphInds)-1
    ind = find(~isnan(kcit(:,1)) & ~isnan(kcit(:,i+1)));
    szCit(i) = numel(ind);
    x1 = kcit(ind,1);
    xn = kcit(ind,i+1);
    diffCit{i} = xn - x1;
    muCit(i) = mean(xn-x1);
    sdCit(i) = std(xn-x1);
    seCit(i) = sdCit(i)/sqrt(szCit(i));
    [h,pvCit(i)] = ttest(x1, xn);
  end
  fprintf('Per/mixture level p-values for citral:  ');
  fprintf('%1.3f ', pvCit);
  fprintf('\n');
  
  mug = sum(cellfun(@sum, diffCit))/sum(cellfun(@numel, diffCit));
  sst = sum(cellfun(@(x) sum((x-mug).^2), diffCit));
  sse = sum(cellfun(@(x) sum((x - mean(x)).^2), diffCit));
  ssr = sum(cellfun(@(x) numel(x)*(mean(x)-mug).^2, diffCit));
  dfe = sum(cellfun(@numel, diffCit)) - numel(diffCit);
  dfr = numel(diffCit) - 1;
  msr = ssr/dfr;
  mse = sse/dfe;
  fstat = msr/mse;
  pvCitOverall = 1 - fcdf(fstat, dfr, dfe);
  
  % Fit a sinusoid
  u = sin(2*pi*(-1.2+(1:10))/9);
  w = [ones(10,1) u(:)]\muCit(:);
  muCite = [ones(10,1) u(:)]*w;
  citRes = muCite(:)-muCit(:);
  r2e = 1 - var(citRes)/var(muCit);
  sse = sum(citRes.^2);
  sst = sum((muCit - mean(muCit)).^2);
  ssr = sst - sse;
  dfr = 1;
  dfe = numel(muCit)-1;
  fstat = (ssr/dfr)/(sse/dfe);
  pval = 1 - fcdf(fstat, dfr, dfe);
  fprintf('Sinusoidal fit of scaling trend for citral:  R2: %1.3e, pval: %1.3e\n', r2e, pval);

  u = sin(2*pi*(-1.5+(1:10))/9);
  muOct = fliplr(muOct);
  w = [ones(10,1) u(:)]\muOct(:);
  muOcte = [ones(10,1) u(:)]*w;
  octRes = muOcte(:)-muOct(:);
  r2e = 1 - var(octRes)/var(muOct);
  sse = sum(octRes.^2);
  sst = sum((muOct - mean(muOct)).^2);
  ssr = sst - sse;
  dfr = 1;
  dfe = numel(muOct)-1;
  fstat = (ssr/dfr)/(sse/dfe);
  pval = 1 - fcdf(fstat, dfr, dfe);
  fprintf('Sinusoidal fit of scaling trend for octanol: R2: %1.3e, pval: %1.3e\n', r2e, pval);
end

%% Metric 5: Is there an effect of morph on the scaling weight Part II - compare to 1, the value we'd have (by definition) at the edges 
if (any(whichMetrics==5))
  % Go through the results and grab, for the citral-type responses
  % the k values of the scaling, and plot them.
  resultsFile = 'resultsBestModel3Bayes1Laplace_slag100_sreg100_maxLag3_funcModelPrior.mat'; %resultsBestModel3Bayes1Laplace_slag100_sreg100_maxLag5_constPrior3.mat';
  Output = LOADF(resultsFile);
  
  morphInds = 1:11;

  R = reshape(Output.Results, 168, []);
  R = R(:,morphInds);

  dataForFigure = LOADF('dataForFigure');
  
  Icit = dataForFigure.I{2};
  Ioct = dataForFigure.I{3};
  
  koct = nan(numel(Ioct), numel(morphInds));
  kcit = nan(numel(Icit), numel(morphInds));;

  R = reshape(Output.Results,168,[]);
  R = R(:,morphInds);
  
  for i = 1:numel(Ioct)
    for j = 1:size(R,2)
      switch(R{Ioct(i),j}.bestModel)
       case {2, 9} % Unit Octanol 
        koct(i,j) = 1;
       case {5,12} % Scaled Octanol
        koct(i,j) = R{Ioct(i),j}.w(2);
      end
    end  
  end
  
  muOct = zeros(1,numel(morphInds));
  sdOct = zeros(1,numel(morphInds));
  seOct = zeros(1,numel(morphInds));
  pvOct = zeros(1,numel(morphInds));
  szOct = zeros(1,numel(morphInds));
  for i = 1:numel(morphInds)
    ind = find(~isnan(koct(:,i)));
    szOct(i) = numel(ind);
    xn = koct(ind,i);
    
    muOct(i) = mean(xn-1);
    sdOct(i) = std(xn-1);
    seOct(i) = sdOct(i)/sqrt(szOct(i));
    [h,pvOct(i)] = ttest(xn);
  end
  
  
  for i = 1:numel(Icit)
    for j = 1:size(R,2)
      switch(R{Icit(i),j}.bestModel)
       case {3, 10} % Unit Citral 
        kcit(i,j) = 1;
       case {6, 13} % Scaled Citral
        kcit(i,j) = R{Icit(i),j}.w(2);
      end
    end
  end
   
  muCit = zeros(1,numel(morphInds));
  sdCit = zeros(1,numel(morphInds));
  seCit = zeros(1,numel(morphInds));
  pvCit = zeros(1,numel(morphInds));
  szCit = zeros(1,numel(morphInds));
  for i = 1:numel(morphInds)
    ind = find(~isnan(kcit(:,i)));
    szCit(i) = numel(ind);
    xn = kcit(ind,i);
    
    muCit(i) = mean(xn-1);
    sdCit(i) = std(xn-1);
    seCit(i) = sdCit(i)/sqrt(szCit(i));
    [h,pvCit(i)] = ttest(xn);
  end
end

%% Metric 6: Compute fraction fit stats at 100:100, because we can
%% compare it to what we find for the equivalent in the complex
%% mixtures setup, odor AC.
if (any(whichMetrics == 6))
  resultsFile = 'resultsBestModel3Bayes1Laplace_slag100_sreg100_maxLag3_funcModelPrior.mat'; %resultsBestModel3Bayes1Laplace_slag100_sreg100_maxLag5_constPrior3.mat';
  Output = LOADF(resultsFile);
  [r2,e2,snr,signalPower, noisePower] = ComputeLinearityIndicesFromBayes1LaplaceResults2(Output.Results, pnCbotDataFile);
  Results = reshape(Output.Results,168,[]);  
  ind100to100 = size(Results,2)-1;
  Results100to100 = Results(:,ind100to100);
  snr100to100 = snr(:,ind100to100);
  snrDb =10*log10(snr100to100);
  nonConst = cellfun(@(x) ~isempty(x.indInputs), Results100to100);
  fracFit = mean(nonConst);
  fracAbove3dB = mean(snrDb>3);
  indAbove3dB = find(snrDb>3);
  fracFitAbove3dB = mean(nonConst(indAbove3dB));
  fprintf('Stats for cit100:oct100, equivalent to odorAC:\n');
  fprintf('Fraction fit: %1.3f\n', fracFit);
  fprintf('Fraction above 3dB: %1.3f\n', fracAbove3dB);
  fprintf('Fraction fit above 3dB: %1.3f\n', fracFitAbove3dB);  
end

%% Metric 7: Fraction of responses to cit100:oct100 that are unit, and lagged
if (any(whichMetrics == 7))
  disp('Fraction of responses to cit100:oct100 that are unit, and lagged');
  resultsFile = 'resultsBestModel3Bayes1Laplace_slag100_sreg100_maxLag3_funcModelPrior.mat'; %resultsBestModel3Bayes1Laplace_slag100_sreg100_maxLag5_constPrior3.mat';
  Output = LOADF(resultsFile);
  Results = Output.Results;
  Results = reshape(Results,168,16);
  Results100to100 = Results(:,end-1);
  bm = cellfun(@(x) x.bestModel, Results100to100);
  bm = bm(:);
  citOrOctFlag = sum(bsxfun(@eq, bm, [2 3 5 6 9 10 12 13]),2)>0;
  unitFlag = sum(bsxfun(@eq, bm, [2 3 9 10]),2)>0;
  lagFlag  = sum(bsxfun(@eq, bm, [9 10 12 13]),2)>0;
  fprintf('CIT100 : OCT100\n');
  fprintf('Fraction of responses cit or oct: %1.3e\n', mean(citOrOctFlag));
  fprintf('Fraction of cit or oct responses with unit fit: %1.3e\n', mean(unitFlag(citOrOctFlag)));
  fprintf('Fraction of cit or oct responses with lag fit: %1.3e\n',  mean(lagFlag(citOrOctFlag)));
end

%% Metric 8: Mean +/sem of scaling weights in preferred cells for preferred responses
if (any(whichMetrics==8))
  disp('Mean +/sem of scaling weights in preferred cells for preferred responses');
  resultsFile = 'resultsBestModel3Bayes1Laplace_slag100_sreg100_maxLag3_funcModelPrior.mat'; %resultsBestModel3Bayes1Laplace_slag100_sreg100_maxLag5_constPrior3.mat';
  Output = LOADF(resultsFile);
  Results = reshape(Output.Results,168,[]);
  Results = Results(:,1:11); % Take the mixture morph results
  dataForFigure = LOADF('dataForFigure');
  indCit = dataForFigure.I{2};
  indOct = dataForFigure.I{3};
  kcit = [];
  koct = [];
  bm = cellfun(@(x) x.bestModel, Results);
  for i = 1:numel(indCit)
    prefFlag = bsxfun(@eq, bm(indCit(i),:), [3 6 10 13]');
    indPref = find(sum(prefFlag,1));
    for j = 1:numel(indPref)
      switch(bm(indCit(i), indPref(j)))
       case {3, 10}
        kcit(numel(kcit)+1) = 1;
       case {6, 13}
        kcit(numel(kcit)+1) = Results{indCit(i), indPref(j)}.w(end);
       otherwise
        error('Expected citral ind prefs to be [3 6 10 or 13].');
      end
    end
  end

  for i = 1:numel(indOct)
    prefFlag = bsxfun(@eq, bm(indOct(i),:), [2 5 9 12]');
    indPref = find(sum(prefFlag,1));
    for j = 1:numel(indPref)
      switch(bm(indOct(i), indPref(j)))
       case {2, 5}
        koct(numel(koct)+1) = 1;
       case {9, 12}
        koct(numel(koct)+1) = Results{indOct(i), indPref(j)}.w(end);
       otherwise
        error('Expected octanol ind prefs to be [2 5 9 or 12].');
      end
    end
  end

  MUSDSE = @(x) [mean(x) std(x) std(x)/numel(x)];
  musdseCit = MUSDSE(kcit);
  musdseOct = MUSDSE(koct);
  musdseBoth = MUSDSE([kcit(:);koct(:)]);
  
  fprintf('Mean +/- SEM, STD for Preferred Responses to Cit: %1.3f +/- %1.3e, %1.3f\n', musdseCit(1), musdseCit(3), musdseCit(2));
  fprintf('Mean +/- SEM, STD for Preferred Responses to Oct: %1.3f +/- %1.3e, %1.3f\n', musdseOct(1), musdseOct(3), musdseOct(2));
  fprintf('Mean +/- SEM, STD for Preferred Responses to Both: %1.3f +/- %1.3e, %1.3f\n', musdseBoth(1), musdseBoth(3), musdseBoth(2));
end

%% Metric 10: What fraction of models for the morph data  were fit with non-constant fits
if (any(whichMetrics==10))
  disp('Fraction of binary mixture morph responses fit with non-constant fits.');
  resultsFile = 'resultsBestModel3Bayes1Laplace_slag100_sreg100_maxLag3_funcModelPrior.mat'; %resultsBestModel3Bayes1Laplace_slag100_sreg100_maxLag5_constPrior3.mat';
  Output = LOADF(resultsFile);
  [r2,e2,snr,signalPower, noisePower] = ComputeLinearityIndicesFromBayes1LaplaceResults2(Output.Results, pnCbotDataFile);
  Results = reshape(Output.Results,168,[]);  
  morphInds = 1:11;
  bm = cellfun(@(x) x.bestModel, Results);
  bm = bm(:,morphInds);
  snrDb = 10*log10(snr(:,morphInds));  
  fracNonConstant = mean(Columnize(bm~=1));
  fracNonConstantAbove3dB = mean(Columnize(bm(snrDb>3)~=1));
  fprintf('Fraction of mixture morph responses fit best with non-constant models: %1.2f\n',       fracNonConstant);
  fprintf('Fraction of mixture morph responses > 3dB fit best with non-constant models: %1.2f\n', fracNonConstantAbove3dB);  
end

%% Metric 11: Fraction of high strong single component responses
if (any(whichMetrics==11))
  data = LOADF('snrAngleData');
  th = data.th;
  ct = data.ct;
  fprintf('Fraction of cit/oct responses within pi/8 of corresponding center angle: %1.3e\n', mean([abs(th(ct==1))<pi/8; abs(th(ct==2)-pi/2)<pi/8]));
end