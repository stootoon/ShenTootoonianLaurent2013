function ProcessData()
% function ProcessData()

figDir     = GetDataDirForFigure(3);
thisDir    = GetCurrentDirFromPathString(fileparts(mfilename('fullpath')));
targetDir  = fullfile(figDir, thisDir, 'recomputedData');

t0             = 1;
binSize        = 0.1;
responseWindow = [1 5];
numBs          = 0; % No bootstrap runs

Mfig = ComputeMetricsForBinaryMixtureTrajectoryProjections(t0, binSize, responseWindow, numBs);

outputFile = fullfile(targetDir, 'Mfig.mat');
save(outputFile, 'Mfig');

fprintf('Wrote "%s".\n', outputFile);
