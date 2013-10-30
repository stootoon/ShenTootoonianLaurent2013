function ProcessData()
% function ProcessData()

global targetDir

figDir    = GetDataDirForFigure(6);
currDir   = GetCurrentDirFromPathString(fileparts(mfilename('fullpath')));
targetDir = fullfile(figDir, currDir ,'recomputedData');

M          = ComputeMetricsForPnKcAuc_noduplicates;
outputFile = fullfile(targetDir, 'Mfig_noduplicates'); 
save(outputFile, 'M');
fprintf('Wrote "%s".\n', outputFile);

binSize    = 0.05;
M          = ComputeMetricsForBinnedAuc_noduplicates(0.05); 
outputFile = fullfile(targetDir, 'Mbin50_noduplicates');
save(outputFile, 'M');
fprintf('Wrote "%s".\n', outputFile);
