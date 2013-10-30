function [logp, D] = ComputePosteriorOfLinearModel(x,y,xr,yr)
% [logp, D] = ComputePosteriorOfLinearModel(x,y,xr,yr)
%
% Computes the posterior probability of a linear model for the 2D data
% whose coordinates are in the vectors X, Y. The model assumes that
%
% y[i] ~ N ( m x[i] + b, sigma^2)
%
% Where M and B are parameters and sigma^2 is the variance of the fit
% residuals. The angle of the slope is assumed to be uniformly
% distributed in [-pi/2 to pi/2], giving a Cauchy distributed prior on
% M. B is assumed to be uniformly distributed in the range given in
% YR.
%
% The posterior is computed by integrating the integrand of this model
% over the parameter ranges in XR and YR. XR and YR each contain the
% min and max ranges of their respective variables. The integral
% itself is estimated using Laplace's trick.
%
% The logarithm of the posterior probability is returned in LOGP, and
% the structure D contains information about the maximum likelihood
% fit.
%
% See also: LinearModelPosteriorIntegrand, EstimateIntegralUsingLaplacesTrick.

k = numel(x);
db = diff(yr);

% Take a guess at the location of the peak by using the values for a linear fit
A = [x(:) ones(numel(x),1)];
w = A\y(:);
a0 = w(1);
b0 = w(2);
v0 = sum((y-a0*x-b0).^2)/(k-1);

% Find the peak numerically
upk = fminsearch(@(u) -log(LinearModelPosteriorIntegrand(u(1),u(2),u(3),x,y,db)), [a0 b0 v0]);
apk = upk(1);
bpk = upk(2);
vpk = upk(3);

% Estimate the probability using Laplace's approximation
logp = log(EstimateIntegralUsingLaplacesTrick(@(u) LinearModelPosteriorIntegrand(u(1),u(2),u(3), x,y,db),[apk, bpk, vpk]));

if (nargout>1) % D was asked for
  D = struct;
  [z, lzpk] = LinearModelPosteriorIntegrand(apk,bpk,vpk, x,y,db); % Get the MAP values 
  D.apk = apk;
  D.bpk = bpk;
  D.vpk = vpk;
  D.lzpk = lzpk;
  D.logp = logp;
  D.fval = @(x) D.apk*x+D.bpk;
end
