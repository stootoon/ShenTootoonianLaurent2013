function D = ComputeCorrelationDistance(X1,X2)
% D = ComputeCorrelationDistance(X1,X2)
% 
% X1, X2 are [numCells x numBins x numOdors] matrices.
% 
% D is a [numBins x numOdors x numBins x numOdors] matrix which
% contains the correlation distances (1 - the correlation) between
% cells in the two inputs. D(odor1, bin1, odor2, bin2) is the
% correlation distance between X1(:, odor1, bin1) and 
% X2(:, odor2, bin2).

[numCells, numBins, numOdors] = size(X1);
D = zeros(numBins, numOdors, numBins, numOdors);

count = 1; totalCount = numOdors*numBins*numOdors;
for i1 = 1:numOdors
  for j1 = 1:numBins
    x1 = X1(:, j1, i1);
    for i2 = 1:numOdors
      ProgressDot2(count,totalCount,numOdors,numBins/2);
      for j2 = 1:numBins
        x2 = X2(:, j2, i2);
        d = 1 - corrcoef([x1 x2]);
        D(j1,i1,j2,i2) = d(1,2);
      end
      count = count + 1;
    end
  end
end

