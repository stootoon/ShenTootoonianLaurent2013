function ax = PlotCorrelationDistanceMatrices(D, spacing)
% ax = PlotCorrelationDistanceMatrices(D, spacing)
%
% Given the numBins x numOdors x numBins x numOdors correlation
% distance matrix D, groups the results by odor comparison pair and
% plots them in a numOdors x numOdors of matrices, each matrix showing
% the numBins x numBins correlation distance plot for that odor
% comparison.  SPACING specifies the horizontal and vertical spacing
% between the plots in the grid. AX contains handles to the axes, so
% that ax(1) is a handle to the plot at (1,1), ax(2) to the plot at
% (1,2), etc.
%
% See also: ComputeCorrelationDistance.

[numBins1, numOdors1, numBins2, numOdors2] = size(D);

minRange = 0;
maxRange = 1;
Q = ComputeSubplotPositions(numOdors1, numOdors2,[],0.05,0.05,spacing,0.05,0.05,spacing);
ax = zeros(numOdors1*numOdors2,1);
nax = 1;
for i = 1:numOdors1
  for j = 1:numOdors2
    ax(nax) = subplotp(Q,GetSubplotIndex(i,j,numOdors1,numOdors2));    
    imagesc(squeeze(D(:,i,:,j)),[minRange, maxRange]);
    set(gca,'xticklabel',[],'yticklabel',[]);
    box off;
    nax = nax+1;
  end
end
    