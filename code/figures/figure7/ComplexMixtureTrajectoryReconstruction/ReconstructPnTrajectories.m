function [Urec,Upn] = ReconstructPnTrajectories(V,U,B,b0)
% [Urec,Upn] = ReconstructPnTrajectories(V,U,B,b0,maxConn,numShuffles)
%
% Given the NUMSAMPLES X NUMKCS matrix V, the NUMSAMPLES x NUMPNS
% matrix U, the NUMKCS x NUMPNS connectivity matrix B, and the 1 x
% NUMPNS offsets vector b0, returns the reconstructions of the PN
% trajectories for the specified odors, in the NUMPNS X NUMBINSPERODOR
% X NUMODORS matrix UREC. The PN trajectories themselves are returned
% in UPN, in the same format.

[numSamples,numPns] = size(U);
[numSamples,numKcs] = size(V);
odors = GetOdorsList;
numOdors = numel(odors);
numBins = numSamples/numOdors;

convert2cbo = @(X) permute(reshape(X,numBins,numOdors,[]),[3 1 2]);
convert2bo_c= @(X) reshape(permute(X,[2 3 1]),numSamples,[]);

Urec = bsxfun(@plus,V*B,b0);
Urec = convert2cbo(Urec);
Upn  = convert2cbo(U);
