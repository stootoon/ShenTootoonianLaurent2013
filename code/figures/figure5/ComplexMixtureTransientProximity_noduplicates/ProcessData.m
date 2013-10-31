function ProcessData()
% function ProcessData()

figDir    = GetDataDirForFigure(5);
currDir   = GetCurrentDirFromPathString(fileparts(mfilename('fullpath')));
targetDir = fullfile(figDir, currDir, 'recomputedData');

numBs            = 0;
t0               = 1;
binSize          = 0.1;
firstBinTime     = 1.5;
whichComparisons = 707;
responseWindow   = [1.5 5];

fprintf('Computing metrics for transient proximity...\n'); tic;
M = ComputeMetricsForTransientProximity_noduplicates(t0, binSize,responseWindow, whichComparisons, numBs);
fprintf('Done in %1.1f seconds.\n', toc);

outputFile = fullfile(targetDir, 'Mfig_noduplicates.mat');
save(outputFile, 'M');
