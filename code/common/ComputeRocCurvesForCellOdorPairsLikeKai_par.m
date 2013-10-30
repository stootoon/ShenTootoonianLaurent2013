function Roc = ComputeRocCurvesForCellOdorPairsLikeKai_par(spikeTimes, r0, r1, verbosity, whichCellsCmps)
% Roc = ComputeRocCurvesForCellOdorPairsLikeKai_par(spikeTimes, r0,r1, verbosity, [whichCells whichCmps])
%
% Given the toc matrix of spiketimes performs ROC analysis on the
% cells by setting a cell's response to be the total number of spike
% across all trials in the [r0, r1] response window.
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

C = CountSpikesInSemiClosedTimeWindow(spikeTimes, r0, r1);

C = reshape(C, numTrials, numOdors, numCells);
C = sum(C);
C = reshape(C, numOdors, numCells)';

odorCmps = 'ABCDWXYZ';
numCmps = numel(odorCmps);

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
auc = zeros(n,1);

parfor i = 1:n
  whichCell = whichCellsCmps(i,1);
  whichCmp  = whichCellsCmps(i,2);
  labels = GetLabelsForAucComputationsForOdorComponent(odorCmps(whichCmp));

  indNzLabels = find(labels~=0);
  labels = labels(indNzLabels);
  scores = C(whichCell, indNzLabels);
  
  cells(i) = whichCell;
  cmps(i) = whichCmp;
  Labels{i} = labels;
  Scores{i} = scores;
  if (numel(unique(scores))>1)
    [x{i},y{i},th,auc(i)] = perfcurve(labels, scores, 1);
  else
    x{i} = [];
    y{i} = [];
    auc(i) = NaN;
  end
end

Roc = struct;
for i = 1:n
  Roc(i).cell = cells(i);
  Roc(i).cmp = cmps(i);
  Roc(i).labels = Labels{i};
  Roc(i).scores = Scores{i};
  Roc(i).x = x{i};
  Roc(i).y = y{i};
  Roc(i).auc = auc(i);
end