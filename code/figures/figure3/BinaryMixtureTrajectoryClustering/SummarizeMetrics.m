function SummarizeMetrics(Outputs)
% function SummarizeMetrics(Outputs)
%
% Writes out a summary of each of the metrics in the cell array
% OUTPUTS. Each element of the cell array should consist of a 2-element cell
% array whose first element is structure returned by COMPUTEMETRICS
% and whose second element is the arguments to COMPUTEMETRICS used to
% compute it.
%
% Outputs will typically have been returned by PARCOMPUTEMETRICS.

for i = 1:numel(Outputs)
  M    = Outputs{i}{1};
  args = Outputs{i}{2};
  name = sprintf('t0=%1.3f bs=%1.3f rw=[%1.1f-%1.1f]', args{1},args{2},args{3}(1),args{3}(2));
  summarizeResults(M, i, name);
end

function summarizeResults(R, id, name)
fprintf('\nSummary for "%s":\n', name);

riByOdor = R.riByOdor;
riByConc = R.riByConc;
prc = [2.5 97.5];
prcvOdor = prctile(riByOdor,prc);
prcvConc = prctile(riByConc,prc);

fprintf('%1.3f - %1.3f CI for rand index by odor: %1.3f - %1.3f.\n', prc(1),prc(2),prcvOdor(1),prcvOdor(2));
fprintf('%1.3f - %1.3f CI for rand index by conc: %1.3f - %1.3f.\n', prc(1),prc(2),prcvConc(1),prcvConc(2));

sfigure(id);
set(gcf,'name', name,'color',[1 1 1]);
clf;
Q = ComputeSubplotPosition(0.05,0.01,0.125,0.05);

subplotp(Q,1);

b = linspace(0,1,100); % histogram bins

fodor = hist(riByOdor,b)/sum(riByOdor); % convert frequency to fraction
hodor = stem(b(fodor>0),fodor(fodor>0),'b','MarkerSize',0,'LineWidth',4); % plot the non-zero bins

hold on;
fconc = hist(riByConc,b)/sum(riByConc,2);
hconc = stem(b(fconc>0),fconc(fconc>0),'r','MarkerSize',0,'LineWidth',4);

xlim([0 1.01]);
set(gca,'FontSize',12);
xlabel('RandIndex','FontSize',14);
ylabel('Fraction Observed','FontSize',14);
title({'Histogram of the bootstrap distribution of the Rand Index', 'for clustering by distance vs clustering by odor (blue)', 'and clustering by distance vs clustering by concentration (red)'},'FontSize',18);
grid on;
