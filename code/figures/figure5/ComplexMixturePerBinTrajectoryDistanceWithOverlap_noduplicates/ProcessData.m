function ProcessData()
% function ProcessData()

figDir    = GetDataDirForFigure(5);
currDir   = GetCurrentDirFromPathString(fileparts(mfilename('fullpath')));
targetDir = fullfile(figDir, currDir, 'recomputedData');

numBs          = 0;
t0             = 1.5;
binSize        = 0.1;
responseWindow = [1.5 5];
numPnShuffles  = 5;

M = ComputeMetricsForPerBinTrajectoryDistanceWithOverlap_nodups(t0, binSize, responseWindow, numBs, numPnShuffles);

outputFile = fullfile(targetDir, 'Mfig_noduplicates.mat');
save(outputFile, 'M');