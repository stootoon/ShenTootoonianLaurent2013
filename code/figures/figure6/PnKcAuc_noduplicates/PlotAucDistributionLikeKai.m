function PlotAucDistributionLikeKai(auc)
% PlotAucDistributionLikeKai(auc)
%
% Plots the auc distribution AUC like Kai does for the figures.

x = 0:0.05:1;
h = hist(auc, x);
barWidth = x(2)-x(1);

X = [x - barWidth/2; x - barWidth/2; x+barWidth/2; x+barWidth/2];
Y = [h; 0*h; 0*h; h];

set(patch(X,Y,[1 0 0]),'FaceColor','none','Edgecolor','blue','LineWidth',0.5);

hold on;

line([0.5;0.5],[0;50],'Color','r','LineWidth',0.5);

mu = mean(auc);

indClosest = argmin(abs(x - mu));

plot(mu, h(indClosest)+3,'v','MarkerSize',5,'MarkerEdgeColor',name2rgb('blue'), 'MarkerFaceColor',name2rgb('blue'));

set(gca,'xtick',[],'ytick',[0 25 50],'yticklabel', [],'XColor',[0.75 0.75 0.75], 'YColor',[0.75 0.75 0.75]);
xlim([-barWidth/2 1+barWidth/2]);
ylim([0 50]);
box on;
axis square;
