function M = ComputeMetricsForGlobalTrajectoryDistanceWithOverlap_nodups(t0, binSize, responseWindow, numBs, baselineWindow)
% M = ComputeMetricsForGlobalTrajectoryDistanceWithOverlap_nodups(t0, binSize, responseWindow, numBs, baselineWindow)
%
% Computes the spearman rank correlation between odor distances and
% trajectory distances. The trajectories are formed by temporally
% concatenating the per-bin trajectories in the response
% window. Trajectory distances are computed both using euclidean and
% correlation distance. Odor distances are measured using overlap
% distance, cosine distance, and another distance measure.
%
% The output is returned in the structure M, which contains fields for
% all the input parameters, and the additional fields below:
%
% SP: a NUMTRAJDISTFUNS x NUMODORDISTFUNS X NUMSAMPLES matrix
% containing the spearman rank correlations for each trial sampling,
% for each pair of trajectory in the response window and odor distance
% functions.
% 
% SPSH: a NUMTRAJDISTFUNS x NUMODORDISTFUNS X NUMSAMPLES matrix
% containing the correlations when the odor labels are shuffled. A new
% shuffling is computed for each trial sampling.
%
% SPBL: a NUMTRAJDISTFUNS x NUMODORDISTFUNS X NUMSAMPLES matrix
% containing the spearman rank correlations for each trial sampling,
% for each pair of trajectories in the baseline window and odor
% distance functions.
%
% TRAJDISTS: a NUMTRAJDISTFUNS cell array of the traj distance functions used.
% ODORDISTS: a NUMODORDISTFUNS cell array of the odor distance functions used.

t1 = CheckStartTimeAndSetEndTimeUsingResponseWindow(t0, [baselineWindow(1), responseWindow(2)], binSize);

Parameters = ComputeInputParametersForTrajectoryComputations_noduplicates('PnComplexMixture','startTime',t0,'endTime',t1,'binSize',binSize,'sampleType', numBs);

B = GetOdorNamesAsBinaryVectors;

cosineDist = 'cosine';
overlapDist      = @(x,Y) arrayfun(@(i) 1 - dot(x,Y(i,:))/max([sum(x), sum(Y(i,:))]), 1:size(Y,1));
intersectionDist = @(x,Y) arrayfun(@(i) 1 - dot(x,Y(i,:))/(sum(x)+sum(Y(i,:)) - dot(x,Y(i,:))), 1:size(Y,1));

odorDistFuns = {overlapDist, intersectionDist, cosineDist};
trajDistFuns = {'correlation', 'euclidean'};

numOdorDists = numel(odorDistFuns);
numTrajDists = numel(trajDistFuns);

t = Parameters.startTime+(0:Parameters.numBins-1)*Parameters.binSize;
responseBins = find(t>=responseWindow(1) & t<responseWindow(2));
baselineBins = find(t>=baselineWindow(1) & t<baselineWindow(2));
whichOdors = [4:11 13:44]; % skip *high and Paraffin oil.
numOdors = numel(whichOdors);

Xresp = Parameters.spikeCounts(:,responseBins,whichOdors, :);
Xresp = reshape(Xresp,[],numOdors,Parameters.numTrials); % Temporally concatenate

Xbl = Parameters.spikeCounts(:,baselineBins,whichOdors, :);
Xbl = reshape(Xbl,[],numOdors,Parameters.numTrials); % Temporally concatenate

Isample = ComputeTrialSamplingFromParametersStructure(Parameters);
numSamples = numel(Isample);
  
sp   = zeros(numTrajDists, numOdorDists, numSamples);
spsh = zeros(numTrajDists, numOdorDists, numSamples); 
spbl = zeros(numTrajDists, numOdorDists, numSamples); 

numDistCmps = numOdors*(numOdors-1)/2;

dresp   = zeros(numTrajDists, numOdorDists, numSamples, numDistCmps);
drespsh = dresp;
dbl     = dresp;
dodors  = dresp;

for i = 1:numTrajDists
  trajDistFun = trajDistFuns{i};
  for j = 1:numOdorDists
    odorDistFun = odorDistFuns{j};
    for k = 1:numSamples
      Yresp = squeeze(mean(Xresp(:,:,Isample{k}),3))';
      Ybl   = squeeze(mean(Xbl(:,:,Isample{k}),3))';
      Yrespsh = squeeze(mean(Xresp(:,randperm(numOdors), Isample{k}),3))'; % shuffle the odor labels
      dresp(i,j,k,:)   = pdist(Yresp,trajDistFun);
      drespsh(i,j,k,:) = pdist(Yrespsh, trajDistFun);
      dbl(i,j,k,:)     = pdist(Ybl,trajDistFun);
      dodors(i,j,k,:)  = pdist(B,odorDistFun);
      sp(i,j,k)    = corr(Columnize(dodors(i,j,k,:)),  Columnize(dresp(i,j,k,:)),   'type', 'Spearman');
      spbl(i,j,k)  = corr(Columnize(dodors(i,j,k,:)),  Columnize(dbl(i,j,k,:)),     'type', 'Spearman');
      spsh(i,j,k)  = corr(Columnize(dodors(i,j,k,:)),  Columnize(drespsh(i,j,k,:)), 'type', 'Spearman');
    end
  end
end

M.t0 = t0;
M.binSize = binSize;
M.responseWindow = responseWindow;
M.baselineWindow = baselineWindow;
M.Isample = Isample;
M.dresp = dresp;
M.drespsh = drespsh;
M.dbl = dbl;
M.dodors = dodors;
M.sp   = sp;
M.spsh = spsh;
M.spbl = spbl;
M.odorDistFuns = odorDistFuns;
M.trajDistFuns = trajDistFuns;