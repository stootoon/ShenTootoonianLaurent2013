function ProcessData()
% function ProcessData()

figDir    = GetDataDirForFigure(3);
currDir   = GetCurrentDirFromPathString(fileparts(mfilename('fullpath')));
targetDir = fullfile(figDir, currDir, 'recomputedData');

addpath(fullfile(GetRootDir(),'code','figures','figure3'));

numMpCores = 6;

t0      = 1.500;
binSize = 0.100;
responseWindow = [1.5 3.5]; %[2.1 3.1];
numBs  = 0;
yrDist = [0 2.0];
yrPAF  = [0 1.5];
nmc    = 10000;

if (MatlabPoolWrapper('size')>0)
  MatlabPoolWrapper('close');
end
MatlabPoolWrapper('open', numMpCores);

%% Compute the results using the log concentration ratio coordinates
coords = 'log'; 
startTime = tic;
Mfig     = ComputeMetricsForBinaryMixtureTrajectoryMorphs(t0, binSize, responseWindow, numBs, yrDist, yrPAF, coords, nmc);
fileName = fullfile(targetDir, ['Mfig' coords '.mat']);
save(fileName, 'Mfig');
fprintf('Wrote "%s" in %1.3f seconds.\n', fileName, toc(startTime));

%% Compute the results using the percentage coordinates
coords    = 'pc'; 
startTime = tic;
Mfig      = ComputeMetricsForBinaryMixtureTrajectoryMorphs(t0, binSize, responseWindow, numBs, yrDist, yrPAF, coords, nmc);
fileName  = fullfile(targetDir, ['Mfig' coords '.mat']);
save(fileName, 'Mfig');
fprintf('Wrote "%s" in %1.3f seconds.\n', fileName, toc(startTime));
