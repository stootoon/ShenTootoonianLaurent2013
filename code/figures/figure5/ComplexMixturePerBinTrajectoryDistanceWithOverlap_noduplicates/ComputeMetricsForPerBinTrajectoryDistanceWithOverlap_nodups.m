function M = ComputeMetricsForPerBinTrajectoryDistanceWithOverlap_nodups(t0, binSize, responseWindow, numBs, numPnShuffles)
% M = ComputeMetricsForPerBinTrajectoryDistanceWithOverlap_nodups(t0, binSize, responseWindow, odorDistFun, trajDistFun, numBs)
%
% Computes the spearman rank correlation between odor distances as and
% trajectory distances for each bin in the response window. The odor
% distance functions used are overlap distance, intersection distance
% and cosine distance, while the trajectory distances used are
% correlation distance and euclidean distance.
%
% The output is returned in the structure M, which contains fields for
% all the input parameters, and the additional fields below:
%
% SP: a NUMTRAJDISTFUNS X NUMODORDISTFUNS X NUMSAMPLES x NUMBINS vector containing the spearman rank correlations
% for each pair of distance functions, each trial sampling, and each bin.
% 
% SPSH: a NUMTRAJDISTFUNS X NUMODORDISTFUNS X NUMSAMPLES x NUMBINS
% vector containing the spearman rank correlations for each pair of
% distance functions, each trial sampling, and each bin, except with
% the odor labels shuffled for the trajectory data.
% 

t1 = CheckStartTimeAndSetEndTimeUsingResponseWindow(t0, responseWindow, binSize);

Parameters = ComputeInputParametersForTrajectoryComputations_noduplicates('PnComplexMixture','startTime',t0,'endTime',t1,'binSize',binSize,'sampleType', numBs);

cosineDist = 'cosine';
overlapDist = @(x,Y) arrayfun(@(i) 1 - dot(x,Y(i,:))/max([sum(x), sum(Y(i,:))]), 1:size(Y,1));
intersectionDist = @(x,Y) arrayfun(@(i) 1 - dot(x,Y(i,:))/(sum(x)+sum(Y(i,:)) - dot(x,Y(i,:))), 1:size(Y,1));

odorDistFuns = {overlapDist, intersectionDist, cosineDist};
trajDistFuns = {'correlation', 'euclidean'};

t = Parameters.startTime+(0:Parameters.numBins-1)*Parameters.binSize;

whichBins = find(t>=responseWindow(1) & t<responseWindow(2));

whichOdors = [4:11 13:44]; % skip *high and Paraffin oil.
odorVectors = GetOdorNamesAsBinaryVectors; % B doesn't contain *high or Paraffin oil, so we don't need to subindex it with whichOdors.

Isample = ComputeTrialSamplingFromParametersStructure(Parameters);
X = Parameters.spikeCounts(:,whichBins,whichOdors, :);
[numCells, numBins, numOdors, numTrials] = size(X);
[sp,dr,dodors, drsh, spsh] = computeCorrelations(X, odorVectors, trajDistFuns, odorDistFuns, Isample);

%% Now compute everything again, but shuffle the PN order
spPnSh = repmat(sp,[1 1 1 1 numPnShuffles]);
drPnSh = repmat(dr,[1 1 1 1 1 numPnShuffles]);
pnShuffles = zeros(numCells, numBins, numOdors, numPnShuffles); % PN identities are shuffled once every bin
for i = 1:numPnShuffles
  Xsh = X;
  for j = 1:numBins
    for k = 1:numOdors
      pnShuffles(:,j,k,i) = randperm(numCells);
      Xsh(:,j,k,:) = X(pnShuffles(:,j,k,i),j,k,:);
    end
  end
  [spPnSh(:,:,:,:,i), drPnSh(:,:,:,:,:,i)] = computeCorrelations(Xsh, odorVectors, trajDistFuns, odorDistFuns, Isample);
end

%% Finally, store everything in a structure to return.
M.t0 = t0;
M.binSize = binSize;
M.responseWindow = responseWindow;
M.Isample = Isample;
M.sp = sp;
M.spsh = spsh;
M.dodors = dodors;
M.dr = dr;
M.drsh = drsh;
M.firstBinTime = t(1);
M.odorDistFuns = odorDistFuns;
M.trajDistFuns = trajDistFuns;
M.numPnShuffles = numPnShuffles;
M.pnShuffles = pnShuffles;
M.spPnSh = spPnSh;
M.drPnSh = drPnSh;

function [sp, dr, dodors, drsh, spsh] = computeCorrelations(X, odorVectors, trajDistFuns, odorDistFuns, Isample)

[numCells, numBins, numOdors, numTrials]  = size(X);
numDistCmps = numOdors*(numOdors-1)/2;
numSamples = numel(Isample);
numOdorDists = numel(odorDistFuns);
numTrajDists = numel(trajDistFuns);

sp = zeros(numTrajDists, numOdorDists, numSamples, numBins);
spsh = sp;

dodors = zeros(numOdorDists, numDistCmps);
dr = zeros(numTrajDists, numOdorDists, numSamples, numBins, numDistCmps);
drsh = zeros(numTrajDists, numOdorDists, numSamples, numBins, numDistCmps);

for i = 1:numTrajDists
  trajDistFun = trajDistFuns{i};
  for j = 1:numOdorDists
    odorDistFun = odorDistFuns{j};
    dodors(j,:) = pdist(odorVectors,odorDistFun);
    for k = 1:numSamples
      Y = mean(X(:,:,:,Isample{k}),4);
      for m = 1:numBins
        dr(i,j,k,m,:)   = pdist(squeeze(Y(:,m,:))',trajDistFun);
        drsh(i,j,k,m,:) = pdist(squeeze(Y(:,m,randperm(numOdors)))',trajDistFun); % shuffle the odor labels
        sp(i,j,k,m)     = corr(Columnize(dodors(j,:)),  Columnize(dr(i,j,k,m,:)),   'type', 'Spearman');
        spsh(i,j,k,m)   = corr(Columnize(dodors(j,:)),  Columnize(drsh(i,j,k,m,:)), 'type', 'Spearman');
      end
    end
  end
end
