function Plot3DTrajectories(Y,binsPerTrajectory,colors)
% Plot3DTrajectories(Y,binsPerTrajectory,colors)
%
% Given the 3 x numTotalBins matrix of points Y, plots trajectories as
% sets of adjacent BINSPERTRAJECTORY points, using the colors given in
% the rows COLORS.

if (size(Y,1)~=3)
  error('Expected Y to have 3 rows.');
end

numPoints = size(Y,2);
if (mod(numPoints, binsPerTrajectory))
  error('Expected the BINSPERTRAJECTORY to divide the number of points.');
end

numTraj = numPoints / binsPerTrajectory;
clf;
for i = 1:numTraj
  inds = (1:binsPerTrajectory)+(i-1)*binsPerTrajectory;
  plot3(Y(1,inds),Y(2,inds),Y(3,inds),'o-','Color',colors(i,:));
  hold on;
end
box on;
axis square;
grid on;
