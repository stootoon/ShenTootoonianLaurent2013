function h = PlotRocLikeKai(X,Y)
% h = PlotRocLikeKai(X,Y)
% 
% Given the cell arrays X and Y of FP and TP data for some cells,
% plots the ROC curves and returns handles to the curves drawn.

for i = 1:numel(X)
  h(i) = plot(X{i}, Y{i},'b');
  hold on;
end

plot([0 1],[0 1],'r','LineWidth',1);
xlim([0 1]);
ylim([0 1]);

set(gca, 'xtick',[],'ytick',[0 0.5 1], 'yticklabel',[], 'XColor',[0.75 0.75 0.75], 'YColor',[0.75 0.75 0.75]);
axis square;
