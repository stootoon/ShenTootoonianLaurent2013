function Roc = ComputeBinnedRocCurvesForCellOdorPairsLikeKai_par(spikeTimes, r0, r1, binSize, verbosity, whichCellsCmps)
% Roc = ComputeRocCurvesForCellOdorPairsLikeKai_par(spikeTimes, r0,r1, binSize, verbosity, [whichCells whichCmps])
%
% Given the toc matrix of spiketimes performs ROC analysis on the
% cells by setting a cell's response to be the total number of spike
% across all trials in the [r0, r1] response window, binned into bins
% of size BINSIZE.
%
% Since ROC analysis is performed only for the mixture experiments, this
% function assumes that there are 44 odors and 7 trials per odor.
%
% verbosity: 'verbose' if a progress bar is desired, 'silent' otherwise.
%
% whichCellsCmps: Empty if auc curves for all cell/odor pairs are
% desired, otherwise a two column matrix, the rows of which specify
% the cell/component pairs to compute the auc distribution for.

numTrials = 7;
numOdors  = 44;
numCells  = size(spikeTimes,2)/numTrials/numOdors;

C = CountSpikesInBinsAndAverageAcrossTrials(spikeTimes,1:numTrials,1:numOdors,1:numCells,'startTime',r0,'endTime',r1,'binSize',binSize);
C = C*numTrials; % To get back the raw counts.
% C is in cbo format. Put it in cob format.
C = permute(C,[1 3 2]); 

odorCmps = 'ABCDWXYZ';
numCmps = numel(odorCmps);

numBins = size(C,3);

load(fullfile(GetRootDir('odors'),'odors.mat'));
load(fullfile(GetRootDir('odors'),'odorschar.mat'));

xvals = linspace(0,1,20);

if (isempty(whichCellsCmps))
  whichCells = repmat(1:numCells, 8, 1);
  whichCmps = repmat((1:8)',numCells,1);
  whichCellsCmps = [whichCells(:) whichCmps(:)];
elseif (size(whichCellsCmps,2)~=2)
  error('Expected whichCellsCmps to be a two-column matrix, but it is not.');
end


n = size(whichCellsCmps, 1);
cells = zeros(n,1);
cmps  = cells;
Labels = {};
Scores = {};
x = {};
y = {};
auc = zeros(n,numBins);

for i = 1:n
  if (verbosity)
    ProgressDot(i,n,25);
  end
  whichCell = whichCellsCmps(i,1);
  whichCmp  = whichCellsCmps(i,2);
  labels = GetLabelsForAucComputationsForOdorComponent(odorCmps(whichCmp));

  indNzLabels = find(labels~=0);
  labels = labels(indNzLabels);
  
  cells(i) = whichCell;
  cmps(i) = whichCmp;
  Labels{i} = labels;
  for j = 1:size(C,3) % Loop over the bins
    scores = C(whichCell, indNzLabels, j);
    Scores{i,j} = scores;
    if (numel(unique(scores))>1)
      [x{i,j},y{i,j},th,auc(i,j)] = perfcurve(labels, scores, 1);
    else
      x{i,j} = [];
      y{i,j} = [];
      auc(i,j) = NaN;
    end
  end
end

Roc = struct;
for i = 1:n
  for j = 1:numBins
    Roc(i,j).cell = cells(i);
    Roc(i,j).cmp = cmps(i);
    Roc(i,j).labels = Labels{i};  
    Roc(i,j).scores = Scores{i,j};
    Roc(i,j).x = x{i,j};
    Roc(i,j).y = y{i,j};
    Roc(i,j).auc = auc(i,j);
  end
end