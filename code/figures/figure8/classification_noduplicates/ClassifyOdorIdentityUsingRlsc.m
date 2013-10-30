function [Results, Scores] = ClassifyOdorIdentityUsingRlsc(Xcbot, lambda)
% [Results, Scores] = ClassifyOdorIdentityUsingRlsc(Xcbot, lambda)
%
% Given the NUMCELLS x NUMBINS x NUMODORS x NUMTRIALS matrix X,
% creates a multi-class regularized least squares classifiers (RLSC)
% by training NUMODORS x NUMODORS-1 binary classifiers in the 'all
% vs. all' configuration. The performance as measured by leave-one-out
% validation for each time bin is returned in the NUMODORS x NUMBINS x
% NUMTRIALS matrix RESULTS, and the raw scores are returned in the
% NUMODORS x NUMODORS x NUMBINS x NUMTRIALS matrix SCORES.

[numCells, numBins, numOdors, numTrials] = size(Xcbot);

Wfull = zeros(numOdors,numCells,numOdors,numBins);
Wpart = zeros(numOdors,numCells,numOdors,numBins,numTrials);

Xobtc = permute(Xcbot,[3 2 4 1]);

parfor i = 1:numOdors
  otherOdors = allbut(1:numOdors,i);
  X1btc = squeeze(Xobtc(i,:,:,:));
  WfullSlice = nan(numCells, numOdors, numBins); % Using NANs so that the i = j cases will stand out.
  WpartSlice = nan(numCells, numOdors, numBins, numTrials);
  for j = otherOdors % The i = j case will not get updated, and will remain with NAN weights.
    X2btc = squeeze(Xobtc(j,:,:,:));    
    for k = 1:numBins
      X = [squeeze(X1btc(k,:,:)); squeeze(X2btc(k,:,:))]; % trials x cells
      y = [ones(numTrials,1);-ones(numTrials,1)];
      WfullSlice(:,j,k) = X'*((X*X'+eye(size(X,1))*lambda)\y);
      for m = 1:numTrials
        X = [squeeze(X1btc(k,allbut(1:numTrials,m),:)); squeeze(X2btc(k,:,:))]; % trials x cells
        y = [ones(numTrials-1,1);-ones(numTrials,1)];
        WpartSlice(:,j,k,m) = X'*((X*X'+eye(size(X,1))*lambda)\y);
      end
    end
  end
  Wfull(i,:,:,:)   = WfullSlice;
  Wpart(i,:,:,:,:) = WpartSlice;
end

Wfull = permute(Wfull,[2 3 1 4]); % c o2 o1 bins
Wpart = permute(Wpart,[2 3 1 4 5]); % c o2 o1 bins trials
Xcobt = permute(Xobtc,[4 1 2 3]); 

%% Now classify the missing trials and record the score.
Results = zeros(numOdors, numBins, numTrials);

% Each column of scores is the classifier scores for the corresponding
% odor. Scores(i,i,b,t) is the own-score, Scores(j,i,b,t) are the
% other classifier scores, for bin b, and for trial t.

Scores  = zeros(numOdors, numOdors, numBins, numTrials);
for i = 1:numOdors
  ProgressDot(i,numOdors,numOdors);
  otherOdors = allbut(1:numOdors,i);
  for j = 1:numBins
    for k = 1:numTrials
      u = Xcobt(:,i,j,k);

      % Wpart(:,:,i,j,k) contains the weights learned for ODORI
      % (without the k'th trial) against all trials of ODOR~I, in time
      % bin J. The dot product > 0 means that this trial is ODORI.

      vpart = squeeze(u' * Wpart(:,:,i,j,k));
      indz = find(vpart == 0);
      vpart(indz) = rand(size(indz)) - 0.5;
      vpart = (vpart>0); % 1 x 43, true if the j'th classifier said that this is odor I.
      ownScore = sum(vpart); % 1 x 1: 
      Scores(i,i,j,k) = ownScore;

      % Wfull(:,:,odor,bin) are the weight learned when using all
      % trials, for classifying ODOR against all others. If the dot
      % product with the first column is > 0, it means that ODOR has
      % been detected.
      WfullSub = Wfull(:,:,:,j);
      
      % We have to replace the learned weights for the subset for
      % which this odor was in the negative class with the (negative)
      % of the weights learned when this trial was missing. Since odor
      % i was the positive class (o1) in these conditions, we have to
      % swap the order of 'otherOdors' and 'i' below.
      WfullSub(:,i,otherOdors) = -Wpart(:,otherOdors,i,j,k);

      % At this point, WfullSub(:,o2,o1) is the weights learned
      % training odor o1 (pos) against odor o2. Those cases in which
      % the present odor showed up as o2 have the partial instead of
      % the full weights.
      
      Bfull = reshape(WfullSub,numCells,[]);
      vfull = reshape(u'*Bfull,numOdors,numOdors); 

      % vfull(o2,o1) is the result from the classifier trained on o1
      % vs o2. The diagonal terms are nans because we didn't learn an
      % odor against itself.
      
      % Randomize the cases that were ties.
      indz = find(vfull == 0);
      vfull(indz) = rand(size(indz)) - 0.5;
      % Find the scores in favour of this odor.
      vfull = (vfull>0); % This will also take care of the nans.
      
      % Get the scores for all classifiers for which this odor was not
      % in the positive class.
      
      otherScores = sum(vfull(:,otherOdors),1); % 1 x numOdors - 1
      Scores(otherOdors,i,j,k) = otherScores;
      
      if (ownScore > max(otherScores))
        Results(i,j,k) = 1;
      elseif (ownScore < max(otherScores))
        Results(i,j,k) = 0;
      else % If it was a tie (and ownscore was non-zero), pick the winner randomly.
        Results(i,j,k) = (ownScore>0)*(randi(sum(ownScore==otherScores))==1);
      end
      
    end
  end
end
  
