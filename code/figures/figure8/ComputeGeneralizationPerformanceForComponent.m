function Results = ComputeGeneralizationPerformanceForComponent(X, whichCmp, numReps, lambda, verbose) 
% Results = ComputeGeneralizationPerformanceForComponent(X, whichCmp, numReps, lambda, verbose) 
%
% SYNPOSIS:
%
%  Given a 4-D matrix of spike counts X whose dimensions correspond to
%  cells, time bins, odors, and trials, computes the generalization
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
%  performed by looping over each of the used odors, removing it, and
%  training the classifier on all remaining trials. The classifier is
%  then validated on all trials of the odor left out. The fraction of
%  trials correctly classified is recorded and stored in RESULTS.
%
%
% INPUTS:
%
%  X: A [cells x time bins x odors x trials] matrix of spike times,
%  t0ypically computed by COUNTSPIKESINBINSANDAVERAGEACROSSTRIALS.
%
%  WHICHCMP: One of 'A','B','C','D','W','X','Y','Z', specifying the
%  component to learn generalization for.
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
%  RESULTS: A [numOdors x numBins x numReps x 2]
%  matrix. RESULTS(I,J,K,1) contains the fraction of validation trials
%  that the classifier correctly classified, when the I'th ingroup
%  odor was removed during training, in the J'th time bin, and the
%  K'th repetition. RESULTS(I,J,K,2) is the same, but for the I'th
%  outgroup odor.
%
% 
% EXAMPLE: Compute generalization performance on A vs ~A, using the PN data.
% 
% pnSpt = LoadTocSpikeTimes('rawpn');
% X = CountSpikesInBinsAndAverageAcrossTrials(pnSpt,map(@Identity,1:7),1:44,pnsToUse,'startTime',2.1,'endTime',4.1,'binSize',0.1);
% Results = ComputeGeneralizationPerformanceForComponent(X,'A', 50, 1, 1);
%
% % Get the coordinates of the outputs, and the outputs themselves.
% 
% [ao,ab,at,av]= ind2sub(size(Results),find(~isnan(Results)));
% r = Results(~isnan(Results));
%
% % Plot the mean performance on the outgroup odor as a function of time bins.
% plot(accumarray(ab(av==2), r(av==2), [], @mean));
%
%
% See also: COUNTSPIKESINBINSANDAVERAGEACROSSTRIALS, GETBALANCEDINGROUPLABELSFORODORCOMPONENT.

[numCells, numBins, numOdors, numTrials] = size(X);

cmpLabels = GetBalancedIngroupLabelsForOdorComponent(whichCmp);

posOdors = find(cmpLabels==1);
negOdors = find(cmpLabels==-1);

minTrials = (min([numel(posOdors) numel(negOdors)])-1)*numTrials;
maxTrials = (max([numel(posOdors) numel(negOdors)])-1)*numTrials;
numReps = when(minTrials == maxTrials,1,numReps);

usedOdors = [posOdors negOdors];
valences  = [ones(numel(posOdors),1); -ones(numel(negOdors),1)];

Results = nan(numOdors,numBins,numReps,2);

for i = 1:numel(usedOdors)
  
  if verbose
    fprintf('Leaving out odor %d of %d.\n',i,numel(usedOdors));
  end

  % Used for indexing into the Results matrix.
  posInClass = when(valences(i)>0, i, i-numel(posOdors));
  classIndex = when(valences(i)>0, 1, 2);
  
  % The positive and negative classes for training.  We need to
  % separate them so that we can sample a fixed number of trials for
  % each of them for training.
  if (i<=numel(posOdors))
    Xp = X(:,:,allbut(posOdors,i),:);
    Xn = X(:,:,negOdors,:);
  else
    Xp = X(:,:,posOdors,:);
    Xn = X(:,:,allbut(negOdors,i-numel(posOdors)),:);
  end

  % The validation trials.
  Xv = X(:,:,usedOdors(i),:);
  
  for j = 1:numBins
    if verbose
      ProgressDot(j,numBins,numBins);
    end
    
    Up = reshape(Xp(:,j,:,:),numCells,[])'; % In trials x cells format
    Un = reshape(Xn(:,j,:,:),numCells,[])';
    
    Uv = reshape(Xv(:,j,:,:),numCells,[])';
    % Don't add a 1's column, Kai doesn't.
    % Uv = [ones(size(Uv,1),1) Uv];
    
    for k = 1:numReps      
      Ip = randperm(size(Up,1));
      In = randperm(size(Un,1));
      
      % The training matrix has MINTRIALS positive trials and negative
      % trials.
      Ut = [Up(Ip(1:minTrials),:);Un(In(1:minTrials),:)];
      
      % Don't add a 1's column, Kai doesn't for his categorization computations.
      % Ut = [ones(size(Ut,1),1) Ut];

      yt = [ones(minTrials,1); -ones(minTrials,1)];

      % Solve for the weights
      % w = (Ut'*Ut + lambda*eye(size(Ut,2)))\(Ut'*yt);
      
      % Use Kai's method for computing the weights
      G = (Ut*Ut' + lambda*eye(size(Ut,1)));
      c = G\yt;
      w = Ut'*c;
      
      % Apply to the validation data
      yv = Uv*w;
      
      % If yv is zero, give it a random +/- value, otherwise leave it.
      yv = 2*((yv+~yv*randn)>0)-1;
      
      Results(posInClass,j,k,classIndex) = mean(yv==valences(i));
    end
  end
end

if verbose
  fprintf('Done.\n');
end