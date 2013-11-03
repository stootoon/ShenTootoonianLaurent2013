function Fit = FitAllModelsForObservations2Laplace(Xobs, yobs, varargin)
% Fit = FitAllModelsForObservations2Laplace(Xobs, yobs, varargin)
%
% SYNPOSIS:
%
%   Given the observed [NUMBINS x NUMPREDICTORS X NUMTRIALS] matrix of
%   regressor data XOBS, regressor and [NUMBINS X NUMTRIALS] response
%   data yobs, performs the fits for all the models, and returns the
%   results in the structure FIT. 
%   
%   The fits are performed by first trial averaging Xobs and yobs to
%   yielding the trial averaged component responses X, and the mixture
%   response y. Log posteriors for the following fits of the mixture
%   response as a function of the component responses are then computed:
%   
%   1. CONSTANT: Mixture response ~ a constant         
%
%      y(t) ~ N(b0, s^2)    
%
%   2. UNIT COMPONENT: Mixture response ~ one of the component responses
%      + offset.
%
%      y(t) ~ N(b0 + x1(t), s^2)
%      y(t) ~ N(b0 + x2(t), s^2)
%
%   3. UNIT MIXTURE: Mixture response ~ the sum of the component
%      responses + offset
%
%      y(t) ~ N(b0 + x1(t) + x2(t), s^2)
%   
%   4. SCALED COMPONENT: Mixture response ~ a scaled version of one of
%      the component responses + offset.
%
%      y(t) ~ N(b0 + b x1(t), s^2)
%      y(t) ~ N(b0 + b x2(t), s^2)
%
%   5. SCALED MIXTURE: Mixture response ~ a scaled version of the sum of
%      the component responses + offset.
%
%      y(t) ~ N(b0 + b (x1(t) + x2(t)), s^2)
%
%   6. FREE MIXTURE: Mixture response ~ an unconstrained linear combination 
%      of the component responses + offset.
%   
%      y(t) ~ N(b0 + b1 x1(t) + b2 x2(t)), s^2)
%
%   7+. LAGGED: Same as above, but allowing the component responses to
%       be lagged. For example, the lagged free mixture model would be of
%       the form
%         
%       y(t) ~ N(b0 + b1 x1(t-t1) + b2 x2(t-t2)), s^2)
%
%   Mixture responses are assumed to be iid normally distributed around
%   their model-specifed mean values, with unknown but fixed noise
%   variance.
%
%   Priors for each of the models are Gaussians centered at zero for
%   the regression coefficients and lag values, and inverse Gamma for
%   the noise variance.
%
%   The fit procedure is very similar for all models. The analytic
%   form of the log posterior integrand is specified as an anonymous
%   function. The peak of this function is then found and its Hessian at
%   the peak estimated numerically, and the integral estimated using
%   Laplace's approximation.
%
%   It would be computationally too expensive to compute posteriors over
%   all possible lag values for each model. As an approximation, for a
%   given model we try all possible lagged combinations of the
%   predictors (within a limited lag range) and find the one that
%   produced that best least square fit to the mixture response. We then
%   lag the predictors accordingly and fit as in the unlagged case, but
%   add a term to the computed posterior to account for the prior on the
%   lags.
% 
% REQUIRED INPUTS:
%
%   XOBS: A NUMTIMEBINS x NUMPREDICTORS x NUMTRIALS matrix containing
%   the binned spike counts of the cell in each trial in response to the
%   component odors.
%
%   YOBS: A NUMTIMEBINS x NUMTRIALS matrix containing the binned spike
%   counts of the cell in each trial in response to the mixture odor.
%
% OPTIONAL INPUTS: 
%
%   The following optional inputs (along with their default values)
%   can be provided as name-value pairs:
%   
%   whichBinsToFit[1:size(Xobs,1)] : The bins to use when fitting models.
%   sigmaReg[4]: The standard deviation of the normal prior on the regression coefficients. 
%   sigmaReg[4]: The standard deviation of the normal prior on the lags.
%   variancePriorAlpha[1]: The alpha parameter of the inverse Gamma prior on the variance.
%   variancePriorBeta[1]:  The beta parameter of the inverse Gamma prior on the variance.
%   lagLimit[5]: The maximum lag value (in magnitude) to check.
%
%   logModelPriors[]: The priors to use for each model. If empty, uses a
%   uniform prior across the models. If not empty, must be either: 
%
%     1) a vector of values, one for each of the models, or 
%     2) A function taking a binary vector specifying which of predictors
%     are being used in each model.
%
%   whichModelsToFit[]: An integer vector specifying the models to
%   fit. If empty, fits all models.
%
% OUTPUTS:
% 
%   A FIT structure has the following fields
%
%   IALL:   The set of all available time bins.
%   IFIT:   The set of time bins used in the fits.
%   START:  The start bin of the bins to fit in the set of all bins.
%   UALL:   Xobs, but in [bins x trials x predictors] format.
%   U:      The subset of the rows of Xobs used to fit.
%   XALL:   A BINS x PREDICTORS matrix dervied from Uall by averaging over trials.
%   X:      The subset of the rows of Xall used to fit.
%   ARGS:   A cell array of arguments passed to each of the fit functions.
%   MODELS: A cell array of structures containing the results of each of the fits.
%   SSE:    A vector of SSE's for each model fit
%   R2:     A vector of R2 vlaues for each model (R2 = SSR/SST).
%   E2:     A vector of E2 vlaues for each model (E2 = SSE/SST).
%   RHO:    A vector of correlation coefficients between the response and each of the fits.
%   LOGMODELLIKELIHOODS: A vector of log likelihoods for each of the models.
%   LOGMODELPOSTERIORS:  A vector of log posteriors  for each of the models.
%   BESTMODEL: The index of the model with the highest posterior.
%   BESTNONCONSTANTMODEL: The index of the non-constant model with the highes

p = inputParser;
p.addOptional('whichBinsToFit',1:size(Xobs,1));
p.addOptional('sigmaReg', 4);
p.addOptional('sigmaLag', 4);
p.addOptional('variancePriorAlpha',1);
p.addOptional('variancePriorBeta',1);
p.addOptional('lagLimit',5);
p.addOptional('logModelPriors',[]);
p.addOptional('whichModelsToFit',[]);
p.parse(varargin{:});

sigmaReg = p.Results.sigmaReg;
sigmaLag = p.Results.sigmaLag;
lagLimit = p.Results.lagLimit;
variancePriorAlpha = p.Results.variancePriorAlpha;
variancePriorBeta  = p.Results.variancePriorBeta;
whichModelsToFit   = p.Results.whichModelsToFit;

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

TrainingFunctions ={@TrainConstantModel,   
                    @TrainConstantModel, % Unit models are trained with the constant model
                    @TrainConstantModel, % because they're equivalent to fitting a constant model to 
                    @TrainConstantModel, % the mixture response after subtracting the component response.
                    @TrainScaledModel,
                    @TrainScaledModel,
                    @TrainScaledModel,
                    @TrainFreeModel,
                    @(X, whichInputs) TrainLaggedModel(@TrainConstantModel, start, nfit, X, whichInputs, lagLimit, 'fast'),
                    @(X, whichInputs) TrainLaggedModel(@TrainConstantModel, start, nfit, X, whichInputs, lagLimit, 'fast'),
                    @(X, whichInputs) TrainLaggedModel(@TrainConstantModel, start, nfit, X, whichInputs, lagLimit, 'fast'),
                    @(X, whichInputs) TrainLaggedModel(@TrainScaledModel,   start, nfit, X, whichInputs, lagLimit, 'fast'),
                    @(X, whichInputs) TrainLaggedModel(@TrainScaledModel,   start, nfit, X, whichInputs, lagLimit, 'fast'),
                    @(X, whichInputs) TrainLaggedModel(@TrainScaledModel,   start, nfit, X, whichInputs, lagLimit, 'fast'),                   
                    @(X, whichInputs) TrainLaggedModel(@TrainFreeModel,     start, nfit, X, whichInputs, lagLimit, 'fast')                    
                    };

numModels = length(TrainingFunctions);

% The inputs to use for each of the models
whichInputs = [0 0; % Constant (no inputs)
               1 0; % Unit octanol 
               0 1; % Unit citral
               1 1; % Unit mixture
               1 0; % Scaled octanol
               0 1; % Scaled citral
               1 1; % Scaled mixture
               1 1; % Free mixture
               1 0; % Same as above, but lagged.
               0 1;
               1 1;
               1 0;
               0 1;
               1 1;
               1 1];

if (size(whichInputs,1)~=numModels)
  error('Number of WHICHINPUTS (%d) does not match number of models (%d).', size(whichInputs,1),numModels);
end

if (isempty(p.Results.logModelPriors))
  logModelPriors = zeros(1,numModels);
elseif (isa(p.Results.logModelPriors, 'function_handle'))
  logModelPriors = arrayfun(@(i) p.Results.logModelPriors(whichInputs(i,:)), 1:size(whichInputs,1));
else
  logModelPriors = p.Results.logModelPriors;
end

if (length(logModelPriors)~=numModels)
  error('Number of models (%d) does not match number of priors provided (%d).', numModels, length(p.Results.logModelPriors));
end

Fit.logModelPriors = logModelPriors;

numInputs = sum(whichInputs,2);

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

% Arguments to pass to the training functions.
% The lagged models need all the bins 
Args = {X,X,X,X,X,X,X,X, Xall, Xall, Xall, Xall, Xall, Xall, Xall};
if (length(Args)~=numModels)
  error('Number of ARGS (%d) does not match number of models (%d).', length(Args), numModels);
end

Fit.Args = Args;
Fit.Models = {};
Fit.sse = [];

if (isempty(whichModelsToFit))
  whichModelsToFit = 1:numModels;
end

for i = 1:numel(whichModelsToFit)
  im = whichModelsToFit(i);
  M = TrainingFunctions{im}(Args{im}, whichInputs(im,:));
  M.descr = ModelDescription(M);
  M.ifitBnds = [ifit(1) ifit(end)];
  M.variancePriorAlpha = variancePriorAlpha;
  M.variancePriorBeta  = variancePriorBeta;
  M.sigmaLag = sigmaLag;
  M.sigmaReg = sigmaReg;
  M.modelId = im;
  Fit.Models{im} = M; 
  Fit.sse(im) = M.sse;
  Fit.r2(im)  = M.r2;
  Fit.e2(im)  = M.e2;
  Fit.rho(im) = M.rho;
  Fit.logModelLikelihoods(im) = M.logModelLikelihood;
  Fit.logModelPosteriors(im)  = Fit.logModelLikelihoods(im) + Fit.logModelPriors(im);
end
bm = whichModelsToFit(argmax([Fit.logModelPosteriors(whichModelsToFit)],1));
Fit.bestModel = bm;
if (bm~=1)
  Fit.bestNonConstantModel = bm;
else
  % Look for a non-constant model
  nonConstantModels = setdiff(whichModelsToFit, 1);
  Fit.bestNonConstantModel = nonConstantModels(argmax(Fit.logModelPosteriors(nonConstantModels)));
end
  
function M = TrainConstantModel(X,whichInputs)
y  = X(:,end);
z  = y - X(:,1:end-1)*whichInputs(:);
v0 = var(z);
b0 = mean(z);
p0 = SanitizeInitialPoint([b0 v0],1e-3);

n  = length(y);
lf = @(bv) ... % bv: last value is the variance, first is the constant value.
     -n/2*log(2*pi*bv(2))...
     -1/2/bv(2)*sum((z-bv(1)).^2)...
     -0.5*log(2*pi*sigmaReg^2)...
     -bv(1)^2/2/sigmaReg^2 ...
     -(variancePriorAlpha+1)*log(bv(2))-variancePriorBeta/bv(2); % Inverse gamma prior (without the constant term)

[Q,xpk] = EstimateIntegralUsingLaplacesTrick(lf, [p0], 'logf');
M.w = xpk(1);
M.v = xpk(2);
M.logModelLikelihood = log(Q);
M.test = @(X) X(:,1:end-1)*whichInputs(:) + M.w;
ye = M.test(X);
res = y - ye;
M.sse = sum(res.^2);
M.e2  = M.sse/(length(res)-1)/var(y);
M.r2  = 1 - M.e2;
M.rho = dot(ye-mean(ye), y-mean(y))/norm(y-mean(y))/norm(ye-mean(ye));
M.lags = 0*whichInputs;
M.indInputs = find(whichInputs);
M.type = 'constant';
M.getX = @(X) X(:,1:end-1);
M.gety = @(X) X(:,end);
end

function M = TrainScaledModel(X,whichInputs)
n = size(X,1);
y = X(:,end);
x = X(:,1:end-1)*whichInputs(:);
lf = @(b0b1v) ... % b0b1v: Last value is the variance, first two values are the constant offset and the scaling value.
     -n/2*log(2*pi*b0b1v(3)) ... % Normalization constant for SSE
     -1/2/b0b1v(3)*sum((y-b0b1v(1)-b0b1v(2)*x).^2) ... % SSE
     -log(2*pi*sigmaReg^2) ... % Normalization constant for regression priors. The constant offset contributes -0.5 log(..), as does the scaling factor => -log(...)
     -(b0b1v(1)^2+b0b1v(2)^2)/2/sigmaReg^2 ... % Regression priors
     -(variancePriorAlpha+1)*log(b0b1v(3))-variancePriorBeta/b0b1v(3); % Inverse gamma prior (without the constant term)

w0 = [ones(size(x)) x]\y;
ye = [ones(size(x)) x]*w0;
v0 = var(ye-y);

p0 = SanitizeInitialPoint([w0' v0], 1e-3);

[Q,xpk] = EstimateIntegralUsingLaplacesTrick(lf, p0, 'logf');
M.w = xpk(1:2);
M.v = xpk(3);
M.logModelLikelihood = log(Q);
M.test = @(X) X(:,1:end-1)*whichInputs(:)*M.w(2) + M.w(1);
ye = M.test(X);
res = y - ye;
M.sse = sum(res.^2);
M.e2  = M.sse/(length(res)-1)/var(y);
M.r2  = 1 - M.e2;
M.rho = dot(ye-mean(ye), y-mean(y))/norm(y-mean(y))/norm(ye-mean(ye));
M.lags = 0*whichInputs;
M.indInputs = find(whichInputs);
M.type = 'scaled';
M.getX = @(X) X(:,1:end-1);
M.gety = @(X) X(:,end);
end

function M = TrainFreeModel(X,whichInputs)
n = size(X,1);
r = sum(whichInputs);
y = X(:,end);
U = [ones(n,1) X(:, find(whichInputs))];

lf = @(bv) ... % bv: last value is the integrand, first value is the constant offset, remainder are the other coefficients.
     -n/2*log(2*pi*bv(r+2)) ... % Normalization constant for SSE
     -1/2/bv(r+2)*sum((y-U*bv(1:r+1)').^2) ... % SSE
     -(r+1)/2*log(2*pi*sigmaReg^2) ... % Normalization constant for regression priors
     -sum(bv(1:r+1).^2)/2/sigmaReg^2 ... % Regression priors
     -(variancePriorAlpha+1)*log(bv(r+2))-variancePriorBeta/bv(r+2); % Inverse gamma prior (without the constant term)

w0 = U\y;
ye = U*w0;
v0 = var(ye-y);

p0 = SanitizeInitialPoint([w0' v0],1e-3);

[Q,xpk] = EstimateIntegralUsingLaplacesTrick(lf, p0, 'logf');
M.w = xpk(1:r+1);
M.v = xpk(end);
M.logModelLikelihood = log(Q);
M.test = @(X) [ones(size(X(:,1))) X(:,find(whichInputs))]*M.w(:);
ye = M.test(X);
res = y - ye;
M.sse = sum(res.^2);
M.e2  = M.sse/(length(res)-1)/var(y);
M.r2  = 1 - M.e2;
M.rho = dot(ye-mean(ye), y-mean(y))/norm(y-mean(y))/norm(ye-mean(ye));
M.lags = 0*whichInputs;
M.indInputs = find(whichInputs);
M.type = 'free';
M.getX = @(X) X(:,1:end-1);
M.gety = @(X) X(:,end);
end

function M = TrainLaggedModel(F,start,n, X, whichInputs, lagLimit, mode)
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

lagVals = CartesianProduct(lags, sum(whichInputs));
metric = 0*lagVals(:,1);
if (isequal(lower(mode),'fast'))
  y = X(start:start+n-1,end);
  % Compute the errors for the least squares fits for all lag combinations
  for ithis = 1:size(lagVals,1)
    U = LagAndExtractColumns(X, indInputs, lagVals(ithis,:), start, n);
    w = [ones(n,1) U(:,1:end-1)]\y;
    ye = [ones(n,1) U(:,1:end-1)]*w;
    metric(ithis) = sum((y-ye).^2);
  end
else
  parfor ithis = 1:size(lagVals,1)
    U = LagAndExtractColumns(X, indInputs, lagVals(ithis,:), start, n);
    Ml = F(U, whichInputs);
    metric(ithis) = -Ml.logModelLikelihood;
  end
end

indBestLags = argmin(metric,1);
bestLags = lagVals(indBestLags,:);
U = LagAndExtractColumns(X, indInputs, bestLags, start, n);
M = F(U, whichInputs);
M.lags = bestLags;
testf = M.test;
M.test = @(X) testf(LagAndExtractColumns(X, indInputs, bestLags, start, n));
M.getX = @(X) LagAndExtractColumns(X(:,1:end-1), indInputs, bestLags, start, n);
M.gety = @(X) X(start+(0:n-1),end);

% Add terms to account for the gaussian prior on the lags
M.logModelLikelihood = M.logModelLikelihood -length(bestLags)/2*log(2*pi*sigmaLag^2) - sum(bestLags.^2)/2/sigmaLag^2;
end
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
    sgn = when(w(i)>0,'+','-');
    mag = abs(w(i));  
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
