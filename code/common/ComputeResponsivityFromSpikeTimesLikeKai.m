function [R,Amp,Rel] = ComputeResponsivityFromSpikeTimesLikeKai(X, dims, baselineStart, responseStart, responseEnd, binWidth, ampTh, relTh)
% [R, Amp, Rel] = ComputeResponsivityFromSpikeTimesLikeKai(X, d, baselineStart, responseStart, responseEnd, binWidth, ampTh, relTh)
%
% Given the toc matrix of spiketimes X with dimensions d, returns a
% binary matrix such that R(i,j) = 1 if cell i was responsive to odor
% j. A cell is said to be responsive to an odor both the amplitude and
% reliability conditions are met:
%
% Amplitude: The average response across trials in at least one of the
% response bins must be at least 'ampTh' baseline stds above the
% basesline mean.  
%
% Reliability: At least relTh/7 trials must contain a
% spike in the response window.
%
% Example: Compute the responsivity of KCs using Kai's parameters.
%
% tocSpikeTimes = LoadTocSpikeTimes('rawkc');
% R = ComputeResponsivityFromSpikeTimesLikeKai(tocSpikeTimes,[7, 44, 209], 1.5,2.1,3.1,0.2,1.5,4);

Amp = ComputeResponseAmplitudeFromSpikeTimesFast(X,   dims, baselineStart, responseStart, responseEnd, binWidth);
Rel = ComputeResponseReliabilityFromSpikeTimes(X, dims, relTh,  responseStart, responseEnd);
R = (Amp>=ampTh) & Rel;