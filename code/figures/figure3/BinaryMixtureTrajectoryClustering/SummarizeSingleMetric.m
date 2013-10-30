function SummarizeSingleMetric(M)
% SummarizeSingleMetric(M)
%
% Summarizes a single instance of the metrics for
% BinaryMixtureTrajectoryClustering.

uconcBl = mean(mean(M.riByConcBl,1),3);
uodorBl = mean(mean(M.riByOdorBl,1),3);
uchanceBl = mean(mean(mean(M.riByChanceBl,1),3),4);  

uconcResp = mean(mean(M.riByConcResp,1),3);
uodorResp = mean(mean(M.riByOdorResp,1),3);
uchanceResp = mean(mean(mean(M.riByChanceResp,1),3),4);  

fprintf('Rand Index for traj dist vs CONC, BASELINE window: ');
mu = mean(uconcBl);
se = sem(uconcBl);
fprintf('%1.3f +/- %1.3f.\n', mu, se);

fprintf('Rand Index for traj dist vs ODOR, BASELINE window: ');
mu = mean(uodorBl);
se = sem(uodorBl);
fprintf('%1.3f +/- %1.3f.\n', mu, se);

fprintf('Rand Index chance level, for BASELINE window: ');
mu = mean(uchanceBl);
se = sem(uchanceBl);
fprintf('%1.3f +/- %1.4f.\n', mu, se);

fprintf('Rand Index for traj dist vs CONC, RESPONSE window: ');
mu = mean(uconcResp);
se = sem(uconcResp);
fprintf('%1.3f +/- %1.3f.\n', mu, se);

fprintf('Rand Index for traj dist vs ODOR, RESPONSE window: ');
mu = mean(uodorResp);
se = sem(uodorResp);
fprintf('%1.3f +/- %1.3f.\n', mu, se);

fprintf('Rand Index chance level, for RESPONSE window: ');
mu = mean(uchanceResp);
se = sem(uchanceResp);
fprintf('%1.3f +/- %1.4f.\n', mu, se);

alpha = 0.001;
h = ttest(uodorResp, uconcResp, alpha, 'right');
if (h)
  fprintf('During the response window, the mean of the RI (odor) > RI (conc) at the %1.3f level.\n', alpha);
else
  fprintf('During the response window, the mean of the RI (odor) IS NOT > RI (conc) at the %1.3f level.\n', alpha);
end


