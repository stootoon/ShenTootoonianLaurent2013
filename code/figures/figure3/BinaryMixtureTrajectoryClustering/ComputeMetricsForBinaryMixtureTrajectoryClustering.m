function M = ComputeMetricsForBinaryMixtureTrajectoryClustering(t0, binSize, responseWindow, numBs, baselineWindow, numKmeans, numChance)
% function M = ComputeMetricsForBinaryMixtureTrajectoryClustering(t0, binSize, responseWindow, numBs, baselineWindow, numKmeans, numChance)
%
% SYNOPSIS:
%
% Computes the trajectory clustering metrics, namely, for each
% trials sample, computes the Rand Index for clustering of
% trajectories in the response window by distance vs by odor. A
% PARAMETERS structure is first computed using the specified input
% arguments.
%
% Rand indices are computed both during the specified baseline and
% response windows. A chance level of the index is computed for both
% regimes as well, by assigning the odor/conc labels (it doesn't
% matter which, since both consist of 9 elements oganized as three
% clusters of three) randomly to the trajectories, and comparing the
% result with the clusterings of trajectories by correlation distance.
%
% INPUTS:
%
% T0: The starting bin to get spike data for. Should be less than the
% start of the baseline window.
% 
% BINSIZE: The binsize to bin the spike counts into.
%
% BASELINEWINDOW/RESPONSEWINDOW: Two element vectors specifying the
% start and ends of the baseline and response windows. Assumed to be
% non overlapping, and for baseline to precede response.
%
% NUMBS: Number of bootstrap runs to use. Set to 0 to just use single trials.
% NUMKMEANS: Number of kmeans runs to perform, since kmeans chooses
% starting points randomly.
% NUMCHANCE: Number of chance runs to perform.
%
% OUTPUT:
%
% In the following, NUMSAMPLES is the number of trials samples
% used. If NUMBS = 0, i.e. if single trials were used, then this will
% equal to the number of trials.
%
% The output structure M contains the following fields:
%
% T0,BINSIZE,BASELINEWINDOW, RESPONSEWINDOW: The input arguments used
% to create the Parameters structure. 
%
% IBS: Containing the indices used to form the bootstrap runs, if
% any. Will be empty of single trials were used.
%
% ICONC: a 10 x 3 matrix containing the indices of the concentrations
% used for each of the 10 concentration subsets.
%
% RIBYCONC[RESP/BL]: A 10 x NUMSAMPLES x NUMKMEANS matrix of rand
% indices comparing trajectories are clustered by conc vs by distance,
% during the response or baseline windows. 10 corresponds to NUMCONC,
% the 10 subsets of 3 of the 5 available concentrations.
%
% RIBYODOR[RESP/BL]: A 10 x NUMSAMPLES x NUMKMEANS matrix of rand
% indices comparing trajectories are clustered by odor vs by
% distance. 10 again is NUMCONC, but all the values across this
% dimension should be the same. The shape is maintained in this way
% for symmetry with RIBYCONC[RESP/BL].
% 
% RIBYCHANCE[RESP/BL]: A 10 x NUMSAMPLES x NUMKMEANS X NUMCHANCE
% matrix containing the rand indices computed by assigning odor/conc
% labels at random to the 9 trajectories for the corresponding kmeans
% run and concentration subset, and comparing against clustering by
% correlation distance, during the baseline or response window.
%
% CLSTR[RESP/BL]BYDISTANCE: A 9 x NUMCONC X NUMSAMPLES X NUMKMEANS
% matrix containing the 9 cluster labels computed for the
% corresponding concentration subset, trials subset, and kmeans run,
% during the response or baseline windows.
%
% CLSTRBYCONCENTRATION: A 9 x NUMCONC x NUMSAMPLES matrix containing the 9
% clustering by concentration labels used for each concentration subset and
% each trial subset. Will vary with the NUMCONC index, but not with NUMSAMPLES.
%
% CLSTRBYODOR: A 9 x NUMCONC x NUMSAMPLES matrix containing the 9
% clustering by odor labels used for each concentration subset and
% each trial subset. Is always [1 1 1 2 2 2 3 3 3]', since the
% concentration or trial subsets don't change the labeling by
% odor. Shaped this way for symmetry with CLSTRBYCONCENTRATION.
%
% NOTES:
%
% UPDATE (Jan 4, 2012) I wasn't happy with the way chance levels were
% being computed, as the mean RI for comparing random labelings of two
% sets of 9 elements into 3 clusters of size 3 each. The reason we did
% this originally is because the odors and concs are by defn already
% in 3 clusters of 3, and the trajectories are (mostly), so the idea
% was to see what labeling we'd get if we just assigne odor/conc
% labels at random. But the trajectories are not always in 3 clusters
% of 3, so a more natural solution is to use the actual trajectory
% clustering results, and assign the odor/conc labels to _these_
% randomly, which is what I've done.
%
% Also, since Kmeans picks random starting points each time, I've
% allowed for doing multiple kmeans runs.
%
% The number of Kmeans and chance calculations are passed in via the
% input parameters NUMKMEANS and NUMCHANCE, respectively.
% 
% Finally, added more documentation.
%
% UPDATE (September 8, 2011): I noticed that if we compare clustering
% by concentration, which has 5 clusters, to clustering by distance,
% which we explicitly set to 3, there'll probably be a negative bias
% in the RI, relative to that for comparing clustering by odor, which
% also has 3. Hence the code is changed to instead look at all 10
% subsets of 3 of the 5 concentrations and compute the rand indices.


t1 = CheckStartTimeAndSetEndTimeUsingResponseWindow(t0, [baselineWindow(1), responseWindow(2)], binSize);

P = ComputeInputParametersForTrajectoryComputations_noduplicates('BinaryMixture','startTime', t0, 'endTime', t1, 'binSize', binSize, 'sampleType', numBs);

% 120's not in the list because there's no 120:120 mixture.
concs = [30 60 80 100 140]'; 

octanolSeries = [concs 0*concs];
citrolSeries  = [0*concs concs];
mixtureSeries = [concs concs];

allSeries = [octanolSeries; mixtureSeries; citrolSeries];

mixtureInds = mapu(@(x) GetIndexForBinaryMixtureConcentrationPair(x(1),x(2)), allSeries');

clstrByConcentration = [1:numel(concs) 1:numel(concs) 1:numel(concs)]';

t = P.startTime + (0:P.numBins-1)*P.binSize;

responseBins = find(t>=responseWindow(1) & t<responseWindow(2));
baselineBins = find(t>=baselineWindow(1) & t<baselineWindow(2));

[foo, Iconc] = CartesianProduct([1:5],3,'Unique');
numUnique = arrayfun(@(i) numel(unique(Iconc(i,:))), 1:size(Iconc,1));
Iconc = Iconc(numUnique==3,:);
numIconc = size(Iconc,1);

if (isempty(P.Ibs)) % Don't bootstrap, just take the single trials.
  Isample = arrayfun(@(i) i, 1:P.numTrials, 'UniformOutput', false);
else
  % Bootstraps were used, so append the bootstrap trial samples and those of the jackknife
  Ibs = mat2cell(P.Ibs,size(P.Ibs,1), ones(1, size(P.Ibs,2)));
  Ijk = arrayfun(@(i) [1:i-1 i+1:P.numTrials], 1:P.numTrials, 'UniformOutput', false);
  Isample = [Ibs Ijk];
end
numSamples = numel(Isample);

riByOdorResp = zeros(numIconc,numSamples, numKmeans);
riByConcResp = zeros(numIconc,numSamples, numKmeans);
riByOdorBl = zeros(numIconc,numSamples, numKmeans);
riByConcBl = zeros(numIconc,numSamples, numKmeans);
riByChanceResp = zeros(numIconc,numSamples, numKmeans, numChance);
riByChanceBl   = zeros(numIconc,numSamples, numKmeans, numChance);

clstrRespByDistance = [];
clstrBlByDistance   = [];
clstrByOdor = [];
clstrByConcentration = [];

L0 = [1 1 1 2 2 2 3 3 3]; % Default labeling used for chance levels.
iprog = 1; % Computation progress indicator
for j = 1:numSamples
  for k = 1:numIconc
    ProgressDot2(iprog, numSamples*numIconc, 1, numIconc);
    % Grab a subset of 3 of the 5 concentrations used.
    concsUsed = concs(Iconc(k,:));
    octInds = arrayfun(@(conc) GetIndexForBinaryMixtureConcentrationPair(conc, 0), concsUsed);
    citInds = arrayfun(@(conc) GetIndexForBinaryMixtureConcentrationPair(0, conc), concsUsed);
    mixInds = arrayfun(@(conc) GetIndexForBinaryMixtureConcentrationPair(conc, conc), concsUsed);    
    allInds = [octInds; citInds; mixInds];
    % Compute the bootstrap trajectory by averaging over the specified
    % trials, and then concatenate into one long vector per mixture.
    respTraj = reshape(mean(P.spikeCounts(:,responseBins, allInds, Isample{j}),4), [], numel(allInds));
    blTraj   = reshape(mean(P.spikeCounts(:,baselineBins, allInds, Isample{j}),4), [], numel(allInds));
    % Cluster the trajectories by odor
    clstrByOdor(:,k,j) = [ones(numel(concsUsed),1); 2*ones(numel(concsUsed),1); 3*ones(numel(concsUsed),1)];
    % Cluster the trajectories by concentration
    clstrByConcentration(:,k,j) = [(1:numel(octInds))'; (1:numel(citInds))' ; (1:numel(mixInds))'];
    for m = 1:numKmeans
      % Cluster the trajectories by distance.
      clstrRespByDistance(:,k,j,m) = kmeans(respTraj', 3, 'Distance', 'correlation'); 
      clstrBlByDistance(:,k,j,m)   = kmeans(blTraj',   3, 'Distance', 'correlation'); 
      % Compute the Rand indices
      riByOdorResp(k,j,m) = ComputeRandIndex(clstrByOdor(:,k,j),          clstrRespByDistance(:,k,j,m));
      riByConcResp(k,j,m) = ComputeRandIndex(clstrByConcentration(:,k,j), clstrRespByDistance(:,k,j,m));
      riByOdorBl(k,j,m)   = ComputeRandIndex(clstrByOdor(:,k,j),          clstrBlByDistance(:,k,j,m));
      riByConcBl(k,j,m)   = ComputeRandIndex(clstrByConcentration(:,k,j), clstrBlByDistance(:,k,j,m));
      % Compute the Chance levels by assigning the odor/conc labels to the trajectories randomly
      for ch = 1:1000
        riByChanceResp(k,j,m,ch) = ComputeRandIndex(L0(randperm(9)), clstrRespByDistance(:,k,j,m));
        riByChanceBl(k,j,m,ch)   = ComputeRandIndex(L0(randperm(9)), clstrBlByDistance(:,k,j,m));
      end
    end
    iprog = iprog + 1;
  end
end

M.riByOdorResp = riByOdorResp;
M.riByConcResp = riByConcResp;
M.riByOdorBl = riByOdorBl;
M.riByConcBl = riByConcBl;
M.riByChanceResp = riByChanceResp;
M.riByChanceBl   = riByChanceBl;
M.t0 = t0;
M.binSize = binSize;
M.responseWindow = responseWindow;
M.baselineWindow = baselineWindow;
M.Ibs = P.Ibs;
M.Iconc = Iconc;
M.clstrRespByDistance = clstrRespByDistance;
M.clstrBlByDistance   = clstrBlByDistance;
M.clstrByOdor = clstrByOdor;
M.clstrByConcentration = clstrByConcentration;
