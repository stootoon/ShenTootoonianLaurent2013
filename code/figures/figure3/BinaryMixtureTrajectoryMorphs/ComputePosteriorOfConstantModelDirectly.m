function [logp, D] = ComputePosteriorOfConstantModelDirectly(x,y,xr,yr)
% [logp, D] = ComputePosteriorOfConstantModelDirectly(x,y,xr,yr)
%
% Computes the posterior probability of a constant model of the 2D
% data whose coordinates are in the vectors X and Y. The model assumes
% that the data are distributed according to
%
% y[i] ~ N(Y1, v).
%
% Y1 is a paramter with uniform prior over the range YR, and V is the
% variance with a Jeffery's prior. Laplace's approximation is used to
% estimate the posterior integral.
%
% The logarithm of the probabilty is returned in LOGP,
% and additional information about the maximum likelihood fit useful
% for plotting is in D.

k = numel(x);
db = diff(yr);

bpk = mean(y); % Value of b at the peak
vpk = sum((y-bpk).^2)/(k+2); % Value of v at the peak

h11 = k/vpk; % -H(1,1) at the peak
h22 = (k+2)/2/vpk^2; % -H(2,2) at the peak.
adh = h11*h22; % Abs(det(H)) at the peak (off-diag values are zero)
q   = 2*pi/sqrt(adh); % Volume of gaussian
fpk = vpk^(-k/2-1)*exp(-(k+2)/2); % value at the peak (without the constant terms)
p   = q*fpk; % Laplace's approximation of the probability

lc = -k/2*log(2*pi) - log(db); % log of the constant terms

logp = log(p) + lc;

if (nargout>1)
  D = struct;
  D.bpk = bpk;
  D.vpk = vpk;
  D.fpk = fpk;
  D.H = -[h11 0; 0 h22];
  D.fval = @(x) D.bpk+0*x;
end

  
