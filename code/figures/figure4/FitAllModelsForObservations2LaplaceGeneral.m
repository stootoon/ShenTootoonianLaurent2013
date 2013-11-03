function [Fit,Models] = FitAllModelsForObservations2LaplaceGeneral(Xobs, yobs, varargin)
% Fit = FitAllModelsForObservations2LaplaceGeneral(Xobs, yobs, varargin)
%
% Just like FITALLMODELSFOROBSERVATIONS2LAPLACE, but accomodates an
% arbitrary number of regressors, and allows for the following
% additional optional parameters:
%
% fitMode['all']: Can be set to 'single' to fit a single model in
% which case the optional parameters WHICHINPUTCONFIG and WHICHMODEL
% must specify the model to fit.
%
% whichInputConfig[]: A binary vector of length equal to the maximum
% number of regressors, indicating those to be used in the fit, when
% the fitMode is 'single'.
%
% whichModel[]: An integer from 1-6 indicating the model to fit when
% the fitMode is 'single'. The values 1-3 indicate the constant,
% scaled, and free models respectively, and the values 4-6 indicate
% the same models but allowing lagged regressors.
%
% The optional argument 'whichModelsToFit' from
% FITALLMODELSFOROBSERVATIONS2LAPLACE is not supported.
%
% See also: FITALLMODELSFOROBSERVATIONS2LAPLACE.

p = inputParser;
p.addOptional('whichBinsToFit',1:size(Xobs,1));
p.addOptional('sigmaReg', 4);
p.addOptional('sigmaLag', 4);
p.addOptional('variancePriorAlpha',1);
p.addOptional('variancePriorBeta',1);
p.addOptional('lagLimit',5);
p.addOptional('logModelPriors',[]);
p.addOptional('fitMode','all');
p.addOptional('whichInputConfig',[]);
p.addOptional('whichModel',[]);
p.parse(varargin{:});

sigmaReg = p.Results.sigmaReg;
sigmaLag = p.Results.sigmaLag;
lagLimit = p.Results.lagLimit;
variancePriorAlpha = p.Results.variancePriorAlpha;
variancePriorBeta  = p.Results.variancePriorBeta;

Context.sigmaReg = sigmaReg;
Context.sigmaLag = sigmaLag;
Context.variancePriorAlpha = variancePriorAlpha;
Context.variancePriorBeta  = variancePriorBeta;

fitMode = p.Results.fitMode;
whichInputConfig = p.Results.whichInputConfig;
whichModel       = p.Results.whichModel;

Fit = struct;

Fit.source = mfilename;

iall = (1:size(Xobs,1))'; % All the bins available
ifit = p.Results.whichBinsToFit(:); % The bins to fit
start= find(iall == ifit(1)); % The index within iall at which ifit starts

Fit.iall = iall;
Fit.ifit = ifit;
Fit.start = start;

nall = length(iall);
nfit = length(ifit);

Fit.logModelPriors = [];

Xobs(:,size(Xobs,2)+1,:) = yobs;

Uall = Xobs;
Uall = permute(Uall,[1 3 2]); % bins x trials x odors

U    = Xobs(ifit,:,:);
U    = permute(U,[1 3 2]); % bins x trials x odors;

Xall = squeeze(mean(Uall,2));
X    = squeeze(mean(U,2));

% The following fields are used by plotting scripts to reproduce the fit.
Fit.Uall = Uall;
Fit.U = U;
Fit.Xall = Xall;
Fit.X = X;

numRegs = size(X,2) - 1;

Fit.logModelPriorFunction = @(modelType, whichInputs, lagged) -sum(whichInputs)*log(numel(whichInputs));
numModels = 1;

switch fitMode
 case 'all'
  numInputConfigs = 2^numRegs;
  
  inputConfigs  = zeros(numInputConfigs, numRegs);
  logLikelihoods= -inf(numInputConfigs, 6);
  logPosteriors = -inf(numInputConfigs, 6);
  exitFlags     = -inf(numInputConfigs, 6);
 
  M = TrainConstantModel(X,zeros(1,numRegs), Context);
  logLikelihoods(1, 1) = M.logModelLikelihood;
  logPosteriors(1, 1)  = logLikelihoods(1, 1) + Fit.logModelPriorFunction(1,zeros(1,numRegs),0);
  
  % Prepare the list of input Args
  allInputConfigs = zeros(numInputConfigs,numRegs+1);
  ii = 2;
  for i = 1:numRegs
    inputCols = nchoosek(1:numRegs,i);
    inputRows = [1:size(inputCols,1)]'*ones(1,size(inputCols,2));
    I = accumarray([inputRows(:) inputCols(:)], 1);
    allInputConfigs(ii+(0:size(I,1)-1),:) = [I i*ones(size(I(:,1)))];
    ii = ii + size(I,1);
  end
  Models = cell(size(allInputConfigs,1),6);
  Models{1,1} = M;

  for i = 2:size(allInputConfigs,1)
    ModelsSlice = cell(1,6);
    whichInputs = allInputConfigs(i,1:end-1);
    
    ModelsSlice{1} = TrainConstantModel(X, whichInputs, Context);
    ModelsSlice{2} = TrainScaledModel(X,   whichInputs, Context);
    ModelsSlice{3} = TrainFreeModel(X,     whichInputs, Context);
    
    ModelsSlice{4} = TrainLaggedModel(@TrainConstantModel, start, nfit, Xall, whichInputs, Context, lagLimit, 'fast');
    % Lag values for other models will be the same, so don't recompute them.
    bestLags       = ModelsSlice{4}.lags; 
    
    ModelsSlice{5} = TrainLaggedModel(@TrainScaledModel,   start, nfit, Xall, whichInputs, Context, lagLimit, 'fast','lagsToUse', bestLags);
    ModelsSlice{6} = TrainLaggedModel(@TrainFreeModel,     start, nfit, Xall, whichInputs, Context, lagLimit, 'fast','lagsToUse', bestLags);
    
    inputConfigs(i,   :) = whichInputs;
    ll = cellfun(@(x) x.logModelLikelihood, ModelsSlice);
    ef = cellfun(@(x) x.exitFlag, ModelsSlice);
    logLikelihoods(i, :) = ll;
    logPosteriors(i,  :) = arrayfun(@(j) ll(j) + Fit.logModelPriorFunction(mod(j-1,3)+1, whichInputs, j>=4), 1:6);
    exitFlags(i, :)      = ef;
    Models(i,:)          = ModelsSlice;
  end  

  %disp('Done.');
  Fit.inputConfigs   = allInputConfigs;
  Fit.logLikelihoods = logLikelihoods;
  Fit.logPosteriors  = logPosteriors;
  Fit.exitFlags      = exitFlags;
  [ii,jj] = find(Fit.logPosteriors==max(Fit.logPosteriors(:)),1);
  Fit.bestInputConfig = allInputConfigs(ii,:);
  Fit.bestModel       = Models{ii,jj};
 case 'single'
  if (~isempty(whichInputConfig) && ~isempty(whichModel))
    if (length(whichInputConfig)~=numRegs)
      error('Size of input configuration must match # of elements.');
    end
    switch (whichModel)
     case 1
      M = TrainConstantModel(X, whichInputConfig, Context);
     case 2
      M = TrainScaledModel(X,   whichInputConfig, Context); 
     case 3
      M = TrainFreeModel(X,     whichInputConfig, Context);
     case 4
      M = TrainLaggedModel(@TrainConstantModel, start, nfit, Xall, whichInputConfig, Context, lagLimit, 'fast');
     case 5
      M = TrainLaggedModel(@TrainScaledModel,   start, nfit, Xall, whichInputConfig, Context, lagLimit, 'fast');
     case 6
      M = TrainLaggedModel(@TrainFreeModel,     start, nfit, Xall, whichInputConfig, Context, lagLimit, 'fast');
     otherwise
      error('Unknown model %d.', whichModel);
    end
    Fit.inputConfigs   = whichInputConfig;
    Fit.logLikelihoods = M.logModelLikelihood;
    Fit.logPosteriors  = Fit.logLikelihoods(1) + Fit.logModelPriorFunction(whichModel, whichInputConfig, whichModel>=4);
    Fit.Model = M;
    Fit.descr = ModelDescription(M);
  else
    error('To fit a single model, but ''whichInputConfig'' and ''whichModel'' must be provided.');    
  end
 otherwise
  error('Unknown fitmode "%s".\n', fitMode);  
end
end % function

function M = TrainConstantModel(X, whichInputs, Context)
y = X(:,end);
z = y - X(:,1:end-1)*whichInputs(:);
v0 = var(z);
b0 = mean(z);
p0 = SanitizeInitialPoint([b0 v0],1e-3);

numUsedRegressors   = sum(whichInputs(:));
numUnusedRegressors = numel(whichInputs) - numUsedRegressors;

n  = length(y);
lf = @(bv) ... % bv(1) = the mean level, bv(2) is the variance
     -n/2*log(2*pi*bv(2))... % Normalization constant for the SSE
     -1/2/bv(2)*sum((z-bv(1)).^2)... % The SSE
     -0.5*log(2*pi*Context.sigmaReg^2)... % Constant for the mean level
     -bv(1)^2/2/Context.sigmaReg^2 ... % Prior on the mean level
     -(Context.variancePriorAlpha+1)*log(bv(2))-Context.variancePriorBeta/bv(2);% Inverse gamma prior (without the constant term)

[Q, xpk, fpk, H, exitFlag] = EstimateIntegralUsingLaplacesTrick(lf, [p0], 'logf');
M.w = xpk(1);
M.v = xpk(2); % The estimated variance
M.exitFlag = exitFlag;
M.logModelLikelihood = log(Q);
M.test = @(X) X(:,1:end-1)*whichInputs(:) + M.w;
ye = M.test(X);
res = y - ye;
M.sse = sum(res.^2);
M.e2  = M.sse/(length(res) - 1)/var(y); 
M.r2  = 1 - M.e2;
M.rho = dot(ye-mean(ye), y-mean(y))/norm(y-mean(y))/norm(ye-mean(ye));
M.lags = 0*whichInputs;
M.whichInputs = whichInputs;
M.indInputs = find(whichInputs);
M.X = X;
M.type = 'constant';
M.getX = @(X) X(:,1:end-1);
M.gety = @(X) X(:,end);
M.descr = ModelDescription(M);
end

function M = TrainScaledModel(X,whichInputs, Context)
n = size(X,1);
y = X(:,end);
x = X(:,1:end-1)*whichInputs(:);

numUsedRegressors   = sum(whichInputs(:));
numUnusedRegressors = numel(whichInputs) - numUsedRegressors;

lf = @(b0b1v) ... % b0b1v(1) is the constant offset, 2 is the scaling factor, 3 is the variance
     -n/2*log(2*pi*b0b1v(3)) ... % Normalization constant for SSE
     -1/2/b0b1v(3)*sum((y-b0b1v(1)-b0b1v(2)*x).^2) ... % SSE
     -log(2*pi*Context.sigmaReg^2) ... % Normalization constant for regression priors. The constant offset contributes -0.5 log(..), as does the scaling factor => -log(...)
     -(b0b1v(1)^2+b0b1v(2)^2)/2/Context.sigmaReg^2 ... % Regression priors
     -(Context.variancePriorAlpha+1)*log(b0b1v(3))-Context.variancePriorBeta/b0b1v(3); % Inverse gamma prior (without the constant term)

w0 = [ones(size(x)) x]\y;
ye = [ones(size(x)) x]*w0;
v0 = var(ye-y);

p0 = SanitizeInitialPoint([w0' v0], 1e-3);
[Q, xpk, fpk, H, exitFlag] = EstimateIntegralUsingLaplacesTrick(lf, p0, 'logf');
M.w = xpk(1:2);
M.v = xpk(3);
M.exitFlag = exitFlag;
M.logModelLikelihood = log(Q);
M.test = @(X) X(:,1:end-1)*whichInputs(:)*M.w(2) + M.w(1);
ye = M.test(X);
res = y - ye;
M.sse = sum(res.^2);
M.e2  = M.sse/(length(res) - 1)/var(y);
M.r2  = 1 - M.e2;
M.rho = dot(ye-mean(ye), y-mean(y))/norm(y-mean(y))/norm(ye-mean(ye));
M.lags = 0*whichInputs;
M.whichInputs = whichInputs;
M.indInputs = find(whichInputs);
M.X = X;
M.type = 'scaled';
M.getX = @(X) X(:,1:end-1);
M.gety = @(X) X(:,end);
M.descr = ModelDescription(M);
end

function M = TrainFreeModel(X,whichInputs,Context)
n = size(X,1);
r = sum(whichInputs);
y = X(:,end);
U = [ones(n,1) X(:, find(whichInputs))];

lf = @(bv) ... % bv(1:end-1) = bv(1) is the constant offset, bv(2:end-1) are the other coefficients, bv(end) is the variance
     -n/2*log(2*pi*bv(r+2)) ... % Normalization constant for SSE
     -1/2/bv(r+2)*sum((y-U*bv(1:r+1)').^2) ... % SSE
     -(r+1)/2*log(2*pi*Context.sigmaReg^2) ... % Normalization constant for regression priors
     -sum(bv(1:r+1).^2)/2/Context.sigmaReg^2 ... % Regression priors
     -(Context.variancePriorAlpha+1)*log(bv(r+2))-Context.variancePriorBeta/bv(r+2); % Inverse gamma prior (without the constant term)

w0 = U\y;
ye = U*w0;
v0 = var(ye-y);

p0 = SanitizeInitialPoint([w0' v0],1e-3);

[Q, xpk, fpk, H, exitFlag] = EstimateIntegralUsingLaplacesTrick(lf, p0, 'logf');
M.w = xpk(1:r+1);
M.v = xpk(end);
M.exitFlag = exitFlag;
M.logModelLikelihood = log(Q);
M.test = @(X) [ones(size(X(:,1))) X(:,find(whichInputs))]*M.w(:);
ye = M.test(X);
res = y - ye;
M.sse = sum(res.^2);
M.e2  = M.sse/(length(res)-1)/var(y);
M.r2  = 1 - M.e2;
M.rho = dot(ye-mean(ye), y-mean(y))/norm(y-mean(y))/norm(ye-mean(ye));
M.lags = 0*whichInputs;
M.whichInputs = whichInputs;
M.indInputs = find(whichInputs);
M.type = 'free';
M.X = X;
M.getX = @(X) X(:,1:end-1);
M.gety = @(X) X(:,end);
M.descr = ModelDescription(M);
end

function M = TrainLaggedModel(F,start,n, X, whichInputs, Context, lagLimit, mode, varargin)
p = inputParser;
p.addOptional('lagsToUse', '');
p.parse(varargin{:});
lagsToUse = p.Results.lagsToUse;

if (sum(whichInputs)<1)
  error('Lagged models require at least one predictor variable.');
end
indInputs = find(whichInputs);
if (start+n-1>size(X,1))
  error('Inconsistent inputs to TrainLaggedModel.');
end
offset = start - 1;
minLag = -min(offset, lagLimit);
maxLag = min(size(X,1) - (offset + n), lagLimit);
lags   = minLag:maxLag;

if (isempty(lagsToUse))
  lagVals = CartesianProduct1(lags, sum(whichInputs));
  metric = 0*lagVals(:,1);
  if (isequal(lower(mode),'fast'))
    y = X(start:start+n-1,end);
    for ithis = 1:size(lagVals,1)
      U = LagAndExtractColumns(X, indInputs, lagVals(ithis,:), start, n);
      w = [ones(n,1) U(:,1:end-1)]\y;
      ye = [ones(n,1) U(:,1:end-1)]*w;
      metric(ithis) = sum((y-ye).^2);
    end
  else
    for ithis = 1:size(lagVals,1)
      U = LagAndExtractColumns(X, indInputs, lagVals(ithis,:), start, n);
      Ml = F(U, whichInputs, Context);
      metric(ithis) = -Ml.logModelLikelihood;
    end
  end
  indBestLags = argmin(metric,1);
  bestLags = lagVals(indBestLags,:);
else
  bestLags = lagsToUse;
end

U = LagAndExtractColumns(X, indInputs, bestLags, start, n);
M = F(U, whichInputs, Context);
M.lags = bestLags;
testf = M.test;
M.test = @(X) testf(LagAndExtractColumns(X, indInputs, bestLags, start, n));
M.getX = @(X) LagAndExtractColumns(X(:,1:end-1), indInputs, bestLags, start, n);
M.gety = @(X) X(start+(0:n-1),end);
M.X = X;
% Add terms to account for the gaussian prior on the lags
M.logModelLikelihood = M.logModelLikelihood -length(bestLags)/2*log(2*pi*Context.sigmaLag^2) - sum(bestLags.^2)/2/Context.sigmaLag^2;
M.descr = ModelDescription(M);
end


function U = LagAndExtractColumns(X,whichColumns,lagAmounts, start, n)
U = X(start:start+n-1,:);
for ithis = 1:numel(whichColumns)
  lag = lagAmounts(ithis);
  U(:,whichColumns(ithis)) = X(lag+(start:start+n-1),whichColumns(ithis));
end
end

function str = ModelDescription(M)
w = M.w;
indInputs = M.indInputs;
lags = M.lags;
str = sprintf('y(t) = %1.3f', w(1));
switch(M.type)
 case 'constant'
  varStr = [];
  for i = 1:numel(indInputs)
    lagStr = when(lags(i),sprintf('%+d',lags(i)),'');
    if (isempty(varStr))
      varStr = sprintf('x_%d(t%s)', indInputs(i), lagStr);
    else
      varStr = sprintf('%s + x_%d(t%s)', varStr, indInputs(i), lagStr);
    end
  end
  if (~isempty(varStr))
    str = [str ' + ( ' varStr ' )'];
  end
 case 'scaled'
  varStr = [];
  for i = 1:numel(indInputs)
    lagStr = when(lags(i),sprintf('%+d',lags(i)),'');
    if (isempty(varStr))
      varStr = sprintf('x_%d(t%s)', indInputs(i), lagStr);
    else
      varStr = sprintf('%s + x_%d(t%s)', varStr, indInputs(i), lagStr);
    end
  end
  str = [str when(w(2)>0,' + ',' - ') sprintf('%1.3f', abs(w(2))) ' ( ' varStr ' )'];
 case 'free'
  for i = 1:numel(indInputs)
    sgn = when(w(i+1)>0,'+','-');
    mag = abs(w(i+1));  
    lagStr = when(lags(i),sprintf('%+d',lags(i)),'');
    str = sprintf('%s %s %1.3f x_%d(t%s)', str, sgn, mag, indInputs(i), lagStr);
  end
 otherwise
  error('Unknown model type "%s".\n', M.type);
end
end

function p0 = SanitizeInitialPoint(p,tol)
% p0 = SanitizeInitialPoint(p,tol)
% 
% Sanitizes the initial point used to search for minima in Laplace
% approximation computations, by setting values that are near zero to
% the (slightly) non-zero value tol. 

p0 = p;
indSmall = find(abs(p)<tol);
p0(indSmall) = sign(p(indSmall))*tol;
p0(abs(p0)<tol/2) = tol; % Values which are exactly zero will still be 0 (since sign(0) = 0), so make these positive.

end
