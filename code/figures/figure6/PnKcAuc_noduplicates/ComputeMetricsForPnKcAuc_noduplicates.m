function M = ComputeMetricsForPnKcAuc_noduplicates
% M = ComputeMetricsForPnKcAuc_noduplicates
%
% Computes the metrics for the PN and KC ROC and AUC curves.
% 
% The returned structures has the following fields:
%
% RPN, RKC: NUMCELLS x NUMODORS binary matrices containing the
% repsonsivity of each cell to each odor.
%
% [PN|KC]CELLS: The cells for the responsive cell-odor pairs.
% [PN|KC]CMPS:  The components (range 1:8) of the odors in the cell-odor pairs.
% ROC[PN|KC]: Structure array of ROC information for the cell odor pairs.
pnSpt = LoadTocSpikeTimes('rawpn');
kcSpt = LoadTocSpikeTimes('rawkc');

pnSpt = ConvertSpikeTimesFromSparseToFull(pnSpt);
kcSpt = ConvertSpikeTimesFromSparseToFull(kcSpt);

Rpn = ComputeResponsivityFromSpikeTimesLikeKai(pnSpt,[7 44 174],1.5,2.1,3.1,0.2,1.5,4);
[pnCells, pnCmps] = find(Rpn(:,4:11));

Rkc = ComputeResponsivityFromSpikeTimesLikeKai(kcSpt,[7 44 209],1.5,2.1,3.1,0.2,1.5,4);
[kcCells, kcCmps] = find(Rkc(:,4:11));

RocPn = ComputeRocCurvesForCellOdorPairsLikeKai_par(pnSpt, 2.1, 3.1, 1, [pnCells pnCmps]);
RocKc = ComputeRocCurvesForCellOdorPairsLikeKai_par(kcSpt, 2.1, 3.1, 1, [kcCells kcCmps]);

aucPn = [RocPn.auc];
aucKc = [RocKc.auc];

M = struct;
M.Rpn = Rpn;
M.Rkc = Rkc;
M.pnCells = pnCells;
M.pnCmps  = pnCmps;
M.kcCells = kcCells;
M.kcCmps  = kcCmps;
M.RocPn = RocPn;
M.RocKc = RocKc;
