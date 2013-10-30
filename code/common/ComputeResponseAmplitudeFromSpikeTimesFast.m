function Amplitude = ComputeResponseAmplitudeFromSpikeTimesFast(spt, dims, baselineStart, responseStart, responseEnd, binSize)
% A = ComputeResponseAmplitudeFromSpikeTimesFast(spt, dims, baselineStart, responseStart, responseEnd, binSize)
%
% Given a toc matrix spike times with trial, odor and cell dimensions in dims, computes,
% for each odor and cell, the mean and std of the spike counts across all baseline bins. 
% It the computes the mean spike count per bin across trials in the response window, and
% computes the maximum across all bins. A(i,j) is the difference between this maximum and 
% the baseline mean normalized by the baseline standard deviation, for cell i in response
% to odor j.
%
% A is the amplitude of the response as defined by Kai.
%
% Example:
%
%  kcSpt = LoadTocSpikeTimes('rawkc');
%  Akc = ComputeResponseAmplitudeFromSpikeTimesFast(tocSpikeTimes, [7 44 209], 1.5, 2.1, 3.1, 0.2);

numTrials = dims(1);
numOdors  = dims(2);
numCells  = dims(3);

[spikeCounts, binStarts] = CountSpikesInBinsAndAverageAcrossTrials(spt, ...
                                                  {1,2,3,4,5,6,7}, 1:numOdors, 1:numCells,...
                                                  'startTime', baselineStart, ...
                                                  'endTime', responseEnd, ...
                                                  'binSize', binSize);

spikeCounts = permute(spikeCounts, [4 3 1 2]); % from cbot to tocb

firstResponseBin = find(binStarts>=responseStart,1);

baselineSpikeCounts = spikeCounts(:,:,:,1:firstResponseBin-1);
responseSpikeCounts = spikeCounts(:,:,:,firstResponseBin:end);

meanBaselineSpikeCounts = zeros(numOdors, numCells);
stdBaselineSpikeCounts  = zeros(numOdors, numCells);
maxOfMeanBinnedResponseSpikeCounts = zeros(numOdors, numCells);

baselineSpikeCounts = permute(baselineSpikeCounts, [1 4 2 3]); % t b o c;
baselineSpikeCounts = reshape(baselineSpikeCounts, [], numOdors, numCells);

% Take the mean and the std of all the baseline bins, over all trials.
meanBaselineSpikeCounts = squeeze(mean(baselineSpikeCounts));
stdBaselineSpikeCounts  = squeeze(std(baselineSpikeCounts));

responseSpikeCounts = permute(responseSpikeCounts, [1 4 2 3]); % t b o c
% Compute the mean spike count in each bin across trials.
meanResponseSpikeCountsPerBinAcrossTrials = squeeze(mean(responseSpikeCounts)); % b o c
% Take the maximum of the the mean spike counts across the bins
maxOfMeanBinnedResponseSpikeCounts = squeeze(max(meanResponseSpikeCountsPerBinAcrossTrials, [], 1)); % o c

Amplitude = ((maxOfMeanBinnedResponseSpikeCounts - meanBaselineSpikeCounts)./stdBaselineSpikeCounts)';



