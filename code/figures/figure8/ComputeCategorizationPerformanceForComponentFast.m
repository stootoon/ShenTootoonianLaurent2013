function [Results, OdorIds] = ComputeCategorizationPerformanceForComponentFast(X, whichCmp, numReps, lambda, verbose) 
% [Results, OdorIds] = ComputeCategorizationPerformanceForComponentFast(X, whichCmp, numReps, lambda, verbose) 
%
% SYNPOSIS:
%
%  Same as COMPUTECATEGORIZATIONPERFORMANCEFORCOMPONENT, but sped up
%  by performing the required matrix divisions using two applications
%  of the Sherman-Morrison rule, rather than the naive approach.
%  Additionally, we return the odor and trial ids used during the
%  categorization. What follows is the help from the naive version.
%
%  Given a 4-D matrix of spike counts X whose dimensions correspond to
%  cells, time bins, odors, and trials, computes the categorization
%  performance of a regularized linear classifier when learning
%  WHICHCMP vs ~WHICHCMP. LAMBDA is the parameter of the ridge
%  regularization. To avoid biasing the classifier, equal numbers of
%  trials are used from the WHICHCMP and ~WHICHCMP classes. If one
%  class is larger than the other, NUMREPS repetitions are performed in
%  which training is performed using random subsets of the larger class
%  of size equal to the smaller class. If VERBOSE is non-zero, the
%  progress of training is printed. 
%
%  Given the component WHICHCMP, the odors are labeled using
%  GETBALANCEDINGROUPLABELSFORODORCOMPONENT. Training is then
%  performed by looping over the bins, and for each repetition,
%  selecting a random subset of trials of the positive and negative
%  class, and for each of these trials, leaving one out, training the
%  classifier, and testing it on the trial left out. The classifier's
%  +/- decision is returned in RESULTS.
%
% INPUTS:
%
%  X: A [cells x time bins x odors x trials] matrix of spike times,
%  typically computed by COUNTSPIKESINBINSANDAVERAGEACROSSTRIALS.
%
%  WHICHCMP: One of 'A','B','C','D','W','X','Y','Z', specifying the
%  component to learn categorization for.
%
%  NUMREPS: The number of repetitions to perform if the number of
%  training trials for the two classes are mismatched.
%
%  LAMBDA: The ridge regression regularization parameter.
%
%  VERBOSE: If non-zero, the status of training will be displayed.
%
%
% OUTPUTS:
%
%  RESULTS: A [numBins x numReps x minTrials x 2]
%  matrix. RESULTS(I,J,K,1) is 1 or 0 based on whether the classifiers
%  decision matched or did not match the true value, for the I'th bin,
%  the J'th repetition, with the K'th positive trial left
%  out. R(I,J,K,2) is the same but for the negative trials.
%
%  ODORIDS: A [numBins x numReps x minTrials x 2 x 2] matrix.
%  RESULTS(I,J,K,1,1) is the odor id of the K'th positive odor used in
%  bin I, rep J. RESULTS(I,J,K,2,1) is the corresponding trial
%  id. Changing the last index to 2 returns the information for the
%  negative odor used.
%
% EXAMPLE: Compute categorization performance on A vs ~A, using the PN data.
% 
% pnSpt = LoadTocSpikeTimes('rawpn');
% X = CountSpikesInBinsAndAverageAcrossTrials(pnSpt,map(@Identity,1:7),1:44,pnsToUse,'startTime',2.1,'endTime',4.1,'binSize',0.1);
% Results = ComputeCategorizationPerformanceForComponent(X,'A',10,1,1);
%
% % Get the coordinates of the outputs, and the outputs themselves.
% 
% [ab, ar, at, ac] = ind2sub(size(Results),find(~isnan(Results)));
% r = Results(~isnan(Results));
%
% % Plot the mean performance on the outgroup odor as a function of time bins.
% plot(accumarray(ab(ac==2), r(ac==2), [], @mean));
%
% See also: COUNTSPIKESINBINSANDAVERAGEACROSSTRIALS, GETBALANCEDINGROUPLABELSFORODORCOMPONENT.

[numCells, numBins, numOdors, numTrials] = size(X);

cmpLabels = GetBalancedIngroupLabelsForOdorComponent(whichCmp);

posOdors = find(cmpLabels==1);
negOdors = find(cmpLabels==-1);

numPos = length(posOdors);
numNeg = length(negOdors);

minTrials = min([numel(posOdors) numel(negOdors)])*numTrials;
maxTrials = max([numel(posOdors) numel(negOdors)])*numTrials;
numReps   = when(minTrials == maxTrials,1,numReps);

Xp = X(:,:,posOdors,:);
Xn = X(:,:,negOdors,:);

Results = nan(numBins, numReps, minTrials, 2);
OdorIds = nan(numBins, numReps, minTrials, 2, 2); % Stores the id and trial for the positive and negative valences for each trial.
tic;
for i = 1:numBins
  
  if verbose
    fprintf('Bin %d / %d.\n', i, numBins);
  end

  Up = reshape(Xp(:,i,:,:),numCells,[])'; % This is now in odors_trials x cells
  Un = reshape(Xn(:,i,:,:),numCells,[])';
  
  for j = 1:numReps

    if (verbose)
      ProgressDot(j,numReps, numReps);
    end
    
    Ip = randperm(size(Up,1));
    In = randperm(size(Un,1));
    Upn = [Up(Ip(1:minTrials),:);Un(In(1:minTrials),:)];
    
    % Rows are in odors_trials format, not trials_odors, so be careful!    
    indp = Ip(1:minTrials);
    trialIdPos = floor((indp-1)/numPos)+1;
    odorIdPos  = posOdors(mod(indp-1,numPos)+1);

    indn = In(1:minTrials);
    trialIdNeg = floor((indn-1)/numNeg)+1;
    odorIdNeg  = negOdors(mod(indn-1,numNeg)+1);

    OdorIds(i,j,:,:,1) = [odorIdPos(:) trialIdPos(:)];
    OdorIds(i,j,:,:,2) = [odorIdNeg(:) trialIdNeg(:)];
    
    % Don't add a 1's column, Kai doesn't.
    % Upn = [ones(size(Upn,1),1) Upn]; 
    ypn = [ones(minTrials,1); -ones(minTrials,1)];

    trialInds = 1:minTrials*2;    

    % Instead of performing a fresh matrix inversion for each trial
    % removed, we're going to speed things up by inverting the
    % unperturbed matrix once, and then computing the trial removed
    % inverses via two applications of the Sherman-Morrison rule.  We
    % can do this because if we think of removing a the e.g. j'th row
    % from a matrix X as actually keeping it inplace but setting it to
    % zero (to yield a matrix J), then 
    %
    % JJ' = XX' - ej*rj - rj'*ej' + ej*rj*ej'
    %
    % where ej is the j'th basis (column) vector, rj is the j'th row of XX'.
    %
    % The first term is then just the j'th row of XX', the second term
    % is the j'th column, and the last is the (j,j)'th element, which
    % we add back to avoid subtracting twice (in the implementation we
    % take care of this by setting rj(j) to zero before making the
    % second update.)
    %
    % The matrix we actually want to invert is (JJ' + lambda I), but
    % this doesn't make a difference to the above.
    
    % Compute the overall inverse
    UU = Upn*Upn';
    G = UU + lambda*eye(size(UU,1));
    iG = inv(G);
    % The terms for the first update we can compute all at once.
    RiG = UU*iG;
    RiG = bsxfun(@rdivide, RiG, 1-diag(RiG));

    % Now compute the solutions for each of the 2*minTrials using two
    % applications of Sherman-Morrison
    for k = 1:minTrials*2
      % First update: 
      % G1 = G - UU(k,:)
      %    = G  - ek x rk
      %    = G + (-ek) x rk
      % So
      % iG1 = iG - iG.(-ek) x rk.iG/(1 + rk.iG.(-ek));
      %     = iG + iG(:,k) x rk.iG/(1 - rk.iG(:,k));
      %

      rk = UU(k,:);
      iG1 = iG + iG(:,k) * RiG(k,:);
      
      % Second update:
      % G2 = G1 - UU(:,k) + UU(k,k)
      %    = G1 - rk'*ek'
      %    = G1 + rk'*(-ek');
      % So
      %
      % iG2 = iG1 - iG1.rk' x (-ek').iG1/(1+(-ek').iG1.rk');
      %     = IG1 + iG1.rk' x iG1(k,:)/(1 - iG1(k,:).rk');
      %     = IG1 + [iG1.rk'/(1 - iG1(k,:).rk')] x iG1(k,:);
      
      rk = rk';
      rk(k) = 0;
      iG1k = iG1*rk;
      iG1k = iG1k/(1-iG1k(k));
      iG2 = iG1 + iG1k*iG1(k,:);
      
      % iG2 = (Ut*Ut + lambda I)^-1, but Ut has the k'th row set to
      % zeros, instead of having it removed. So to get the weights, we
      % multiply yt by this, after setting the k'th element of yt to
      % zero. We then grab all but the k'th element of the resulting
      % vector.

      inds = allbut(trialInds, k);
      
      Ut = Upn(inds,:);      
      yt = ypn;
      yt(k) = 0;

      %G = (Ut*Ut'+lambda*eye(size(Ut,1)));
      %c = G\yt;
      c = iG2*yt;      
      w = Ut'*c(inds);
      
      Uv = Upn(k,:);            
      yv = Uv*w;
      
      % If yv~=0, the class is its sign. if yv == 0, set it randomly to +/-.
      yv = 2*((yv + ~yv*randn)>0)-1;      
      
      % Valence-based indexes into the Results matrix.
      posInClass = mod(k-1,minTrials)+1;
      classIndex = 1+(k>minTrials);
      
      Results(i,j,posInClass, classIndex) = (yv == 2*(k<=minTrials)-1);
    end
  
  end

  if verbose
    elapsed = toc;
    eta = elapsed/i*(numBins-i);
    fprintf('~%d seconds remaining.\n', round(eta));
  end

end

if verbose
  fprintf('Done.\n');
end

