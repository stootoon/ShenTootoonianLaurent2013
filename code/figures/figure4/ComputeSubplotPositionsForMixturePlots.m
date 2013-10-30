function Q = ComputeSubplotPositionsForMixturePlots(numPns,numMixtures, horizontalSpacing, maxHeight)
% Q = ComputeSubplotPositionsForMixturePlots(numPns,numMixtures, horizontalSpacing, maxHeight)
%
% Given the vector of mixture responses to plot for the specified
% number of PNs, each response matrix having NUMMIXTURES columns, and
% with the matrices spaced HORIZONTALSPACING apart, computes subplot
% positions that preserve the 1:1 PN to mixture aspect ratio. If the
% initial height is too large and some plots fall outside figure, the
% maximum height is reduced until the plots fit.

n = numel(numPns);
maxPns = max(numPns);
ar = numMixtures/maxPns;
while(1)
  unitWidth = ar*maxHeight;
  innerWidth = (n-1)*(unitWidth+horizontalSpacing); % Width between first and last plot centers
  unitHeights = numPns/maxPns*maxHeight;
  leftMostLeft = 0.5 - innerWidth/2 - unitWidth/2;
  plotLefts = leftMostLeft+(0:n-1)*(unitWidth+horizontalSpacing);
  plotBottoms = 0.5 - unitHeights/2;
  if (any(plotLefts<0))
    maxHeight = maxHeight*0.99;
  else
    break;
  end
end

Q = [plotLefts(:) plotBottoms(:) unitWidth*ones(n,1) unitHeights(:)];



