function Ush = CreateFakeCellsByShufflingOdorResponses(U)
% Ush = CreateFakeCellsByShufflingOdorResponses(U)
%
% Given the CBO or CBOT matrix of cell responses U, creates a set of
% fake cells by shuffling the odor respones, while keeping the trial
% responses constant.

switch(numel(size(U)))
 case 3
 [numCells,numBins,numOdors] = size(U);
 Ush = 0*U;
 for i = 1:numOdors
   Ush(:,:,i) = U(randperm(numCells),:,i);
 end
 case 4
  [numCells,numBins,numOdors,numTrials] = size(U);
  Ush = 0*U;
  for i = 1:numOdors
    Ush(:,:,i,:) = U(randperm(numCells),:,i,:);
  end
 otherwise
  error('U should have 3 or 4 dimensions.');
end


