function [z,lz,s2] = LinearModelPosteriorIntegrand(a,b,v,x,y,db)
% [z,lz] = LinearModelPosteriorIntegrand(a,b,v,x,y,db)
%
% Computes the integrand used for computing the likelihood of a linear
% model to fit the 2D data in the vectors X and Y. The integrand is of the form 
%
% P(D|a,b,v,model)P(a,b|model) = P(D|a,b,v,model)P(a|model)P(b|model)P(v|model)
%
% where a is the slope of the fit, b is the offset, and v is the
% variance. The slope, offset, and variance are assumed to be
% independent, so the prior on the parameters factorizes as
% P(a|model)P(b|model)P(v|model).
%
% The model assumes that the data points are independent and
% distributed according to y[i] ~ N( a x[i] + b, v). The prior on b is
% assumed to be uniform and has width DB. The prior on the slope A is
% taken so that the _angle_ of the slope is uniformly distributed in
% [-pi/2 pi/2]. This yields a prior on the slope of 1/(pi (a^2 + 1)),
% a Cauchy distribution. The prior on the variance v is taken to be a
% Jeffery's prior of 1/v.
%
% A should be a row vector slopes. B should be a single offest
% value. V should be a single variance value. X and Y should be the x-
% and y- coordines of the data points to fit. DB is the length of the
% prior range on B.
%
% Z is the resulting vector of integrand values, and LZ is their natural logarithm.
% 
% This model is suitable for passing to DBLQUAD. Note that if there
% are many data points (>=O(100)), the DBLQUAD result may be inaccurate.

a = a(:)';
k = numel(x);
res = bsxfun(@minus, y(:), x(:)*a+b);
R2 = sum(res.^2,1);
lz = -k/2*log(2*pi)-log(db*pi)-log(1+a.^2)-(k+2)/2*log(v)-R2./(2*v);
z = exp(lz);
