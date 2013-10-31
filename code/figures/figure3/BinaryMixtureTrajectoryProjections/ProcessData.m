function ProcessData()
% function ProcessData()
whichFigure = 3;

figDir     = GetDataDirForFigure(whichFigure);
thisDir    = GetCurrentDirFromPathString(fileparts(mfilename('fullpath')));
targetDir  = fullfile(figDir, thisDir, 'recomputedData');

addpath(GetCodeDirForFigure(whichFigure)); % Add figure3 because it has ComputePafOverallAndPerBin

t0             = 1;
binSize        = 0.1;
responseWindow = [1 5];
numBs          = 0; % No bootstrap runs

Mfig = ComputeMetricsForBinaryMixtureTrajectoryProjections(t0, binSize, responseWindow, numBs);

outputFile = fullfile(targetDir, 'Mfig.mat');
save(outputFile, 'Mfig');

fprintf('Wrote "%s".\n', outputFile);
