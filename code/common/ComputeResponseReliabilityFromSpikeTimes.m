function [R,N] = ComputeResponseReliabilityFromSpikeTimes(spt, dims, rth,t0,t1)
% [R,N] = ComputeResponseReliabilityFromSpikeTimes(spikeTimes, dims, rth, t0, t1)
% 
% [R,N] = ComputeKcReliabilityFromSpikeTimes(spikeTimes, rth, t0, t1)
% 
% Takes a toc matrix of KC spike times SPT with dimensions DIMS and
% returns the [nKcs x nOdors] reliability matrix R. R is a binary
% matrix with r_ij set to 1 if KC i produced spikes in the [t0 t1]
% interval on rth or more trials in response to odor j. N has the same
% format as R but contains the total number of responsive trials. So
% R(i,j) = N(i,j)>=rth
%
% This function assumes that there are 7 trials per odor and 44
% odors.

numTrials = dims(1);
numOdors  = dims(2);
numCells  = dims(3);

nr = CountSpikesInSemiClosedTimeWindow(spt,t0,t1) > 0;

% nr is now an indicator matrix for spikes in [t0 t1] for each trial
nr = reshape(nr,numTrials, numOdors*numCells);
nr = sum(nr,1);
N = reshape(nr, numOdors,numCells)';

% N is now a matrix containing the number of responsive trials in
% the time window.
R = N>=rth;



