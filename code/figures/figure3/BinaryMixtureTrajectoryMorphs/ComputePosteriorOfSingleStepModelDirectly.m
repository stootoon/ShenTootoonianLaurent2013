function [logp, D] = ComputePosteriorOfSingleStepModelDirectly(x,y,xr,yr)
% [logp, D] = ComputePosteriorOfSingleStepModelDirectly(x,y,xr,yr)
%
% Computes the posterior of a Single Step model of the 2D data whose
% coordinates are in the vectors X and Y. The model assumes that the
% data are distributed according to
%
% y[i] ~ N( Y1, sigma^2) : If x[i] <  TH
% y[i] ~ N( Y2, sigma^2) : If x[i] >= TH
% 
% Y1 is the subthreshold level, Y2 is the suprathreshold level, and TH
% is the threshold. All three are parameters. The priors on Y1 and Y2
% is uniform in the range specified by YR, and that on the threshold
% is uniform in the range specified by XR. 
%
% The posterior is computed by integrating the product of the
% likelihood and the prior over the parameter ranges in XR and YR. The
% integral is approximated using Laplace's trick, using analytically
% derived values of the Hessian. The prior on the threshold to is
% limited to be within the range of x values received, so that all
% models are in fact step models.
%
% The logarithm of the posterior probability is returned in LOGP, and
% information about the maximum likelihood fit, useful for plotting,
% is returned in the stucture D.
%
% See also: StepModelPosteriorIntegrand, EstimateIntegralUsingLaplacesTrick.

% The integration takes advantage of the fact that the value of the
% integrand as a function of theta is constant in between datapoints,
% since the set of sub- and super-threshold elements doesn't change.
% Hence we just move from one element to the next, compute the value
% of the integral once, and multiply by the sizes of the x-values
% between the data points.

dx = diff(x); % length of spaces between datapoints.
db = diff(yr);
k  = numel(x);
k1 = 1:numel(x)-1; % # subthreshold elements for each step.
k2 = numel(x)-k1;  % # suprathreshold elements for each step.

% Make space for the variables.
y1 = zeros(1,k-1); % peak values of the first step
y2 = zeros(1,k-1); % peak values of the second step
r1 = zeros(1,k-1); % sum of sub-threshold residuals
y2 = zeros(1,k-1); % sum of supra-threshold residuals.

for i = 1:k-1
  y1(i) = mean(y(1:i));
  y2(i) = mean(y(i+1:end));
  r1(i) = sum((y(1:i)-y1(i)).^2);
  r2(i) = sum((y(i+1:end)-y2(i)).^2);
end

v   = (r1+r2)./(k+2); % Value of v at the peak.
h11 = k1./v; % Diagonal elements of the Hessian. 
h22 = k2./v;
h33 = (k+2)./(2*v)./v;
adh = h11.*h22.*h33; % abs(det(Hessian at the peak)). Off-diagonal elements are zero.
q   = (2*pi)^(3/2)./sqrt(adh); 
fpk = v.^(-k/2-1).*exp(-(k+2)/2);
pl  = q.*fpk; % laplace's approximation
p   = dot(pl,dx); 

% The constant part of the integral. Note that we ignore the passed in
% value of XR, and use the actual values of X that we have. This is
% because we only allow the threshold to lie within our x domain so
% that we don't have degenerate cases of constants etc on the ends.
lc = -k/2*log(2*pi)-2*log(db)-log(x(end)-x(1));

logp = log(p)+lc;

if (nargout>1) % D was asked for
  % Fill D with the values of the fit with the highest likelihood.
  D = struct;
  imax = argmax(pl);
  D.a  = y1(imax);
  D.b  = y2(imax);
  D.v  = v(imax);
  D.th = x(imax+1);
  D.fval = @(x) (x<D.th)*D.a + (x>=D.th)*D.b;%% Make the decision. 
end
