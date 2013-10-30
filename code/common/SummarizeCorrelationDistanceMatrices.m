function Ds = SummarizeCorrelationDistanceMatrices(D, f)
% Ds = SummarizeCorrelationDistanceMatrices(D, f)
%
% Collapses the full correlation distance matrix to a 
% numOdors1 x numOdors2 summary matrix where Ds(i,j) = f(D(:,i,:,j));


[numBins1, numOdors1, numBins2, numOdors2] = size(D);

Ds = zeros(numOdors1, numOdors2);

for i = 1:numOdors1
  for j = 1:numOdors2
    Ds(i,j) = f(squeeze(D(:,i,:,j)));
  end
end