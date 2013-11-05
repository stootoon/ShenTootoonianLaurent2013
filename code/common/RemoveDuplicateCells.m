function [Y, dims, duplicates, removed] = RemoveDuplicateCells(X, dims, corrThresh)
% [Y, dims, duplicates, removed] = RemoveDuplicateCells(X, dims, corrThresh)
%
% Given the TOC matrix X with dimensions DIMS, keeps only 1 example
% whenever the correlation is >= corrThresh (defaults to 1 - 1e-3). A
% sparse toc matrix Y is returned with updated DIMS after removing the
% duplicates. The duplicates themselves are returned as a cell array
% DUPLICATES, and the indices of the cells removed are in REMOVED.

if (nargin==2)
  corrThresh = 1 - 1e-3;
end

maxSpikes = size(X,1);
X = ConvertSpikeTimesFromSparseToFull(X);
X(isinf(X)) = sqrt(-1);

numCells = dims(end);
X = reshape(X,[],numCells);

C = corrcoef(real(X));

duplicates = {};

for i = 1:numCells
  indOthers = i+1:numCells;
  indDup = find(C(i,indOthers)>=corrThresh);
  if (~isempty(indDup))
    duplicates{length(duplicates)+1} = [i indOthers(indDup)];
  end
end

% Check that the duplicates are disjoint
for i = 1:numel(duplicates)
  I = duplicates{i};
  for j = i+1:numel(duplicates)
    J = duplicates{j};
    if (~isempty(intersect(I,J)))
      error('Duplicates are not disjoint.');
    end
  end
end

Y = reshape(X,[],dims(1),dims(2),dims(3));

% Remove the duplicates
removed = [];
for i = 1:numel(duplicates)
  removed = [removed duplicates{i}(2:end)];
end

if (~isempty(removed))
  fprintf('Removing duplicates PNs. IDs of PNs removed: ');
  fprintf('%d ', removed);
  fprintf('\n');
else
  fprintf('No duplicates found.\n');
end

indCells = 1:numCells;
indRemaining = setdiff(indCells, removed);

Y = Y(:,:,:,indRemaining);
Y = reshape(Y,maxSpikes,[]);
Y = sparse(real(Y));
dims(end) = numel(indCells);
