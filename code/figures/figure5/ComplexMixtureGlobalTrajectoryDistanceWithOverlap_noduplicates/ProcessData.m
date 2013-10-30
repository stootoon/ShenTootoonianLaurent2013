function ProcessData()
% function ProcessData()

numBs          = 0;
t0             = 0.5;
binSize        = 0.1;
responseWindow = [2.1 3.1];
baselineWindow = [0.9 1.9];
numBs          = 0;

M = ComputeMetricsForGlobalTrajectoryDistanceWithOverlap_nodups(t0, binSize, responseWindow, numBs, baselineWindow);

figDir     = GetDataDirForFigure(5);
currDir    = GetCurrentDirFromPathString(fileparts(mfilename('fullpath')));
targetDir  = fullfile(figDir, currDir, 'recomputedData');
outputFile = fullfile(targetDir, 'Mfig_noduplicates.mat');
save(outputFile, 'M');