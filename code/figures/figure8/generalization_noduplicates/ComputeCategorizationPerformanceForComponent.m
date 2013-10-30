function Results = ComputeCategorizationPerformanceForComponent(X, whichCmp, numReps, lambda, verbose) 
% Results = ComputeCategorizationPerformanceForComponent(X, whichCmp, numReps, lambda, verbose) 
%
% SYNPOSIS:
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

minTrials = min([numel(posOdors) numel(negOdors)])*numTrials;
maxTrials = max([numel(posOdors) numel(negOdors)])*numTrials;
numReps   = when(minTrials == maxTrials,1,numReps);

Xp = X(:,:,posOdors,:);
Xn = X(:,:,negOdors,:);

Results = nan(numBins, numReps, minTrials, 2);
tic;
for i = 1:numBins
  
  if verbose
    fprintf('Bin %d / %d.\n', i, numBins);
  end

  Up = reshape(Xp(:,i,:,:),numCells,[])';
  Un = reshape(Xn(:,i,:,:),numCells,[])';
  
  for j = 1:numReps

    if (verbose)
      ProgressDot(j,numReps, numReps);
    end
    
    Ip = randperm(size(Up,1));
    In = randperm(size(Un,1));
    Upn = [Up(Ip(1:minTrials),:);Un(In(1:minTrials),:)];
    
    % Don't add a 1's column, Kai doesn't.
    % Upn = [ones(size(Upn,1),1) Upn]; 
    
    ypn = [ones(minTrials,1); -ones(minTrials,1)];

    trialInds = 1:minTrials*2;    

    for k = 1:minTrials*2
      Ut = Upn(allbut(trialInds, k),:);      
      yt = ypn(allbut(trialInds, k),:);

      % w = (Ut'*Ut+lambda*eye(size(Ut,2)))\(Ut'*yt);

      % Use Kai's method for computing the weights.
      G = (Ut*Ut'+lambda*eye(size(Ut,1)));
      c = G\yt;
      w = Ut'*c;
      
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

