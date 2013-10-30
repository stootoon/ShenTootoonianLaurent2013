function ProcessData()
% This script calls COMPUTEMETRICS with different input arguments and
% stores the results in the cell array M, whose elements are each 2
% element cells whose first element is the metrics structure computed
% by COMPUTEMETRICS for a given set of arguments, and the second
% element is the arguments used. The results are typically passed to
% SUMMARIZEMETRICS.

% Format of the arguments below:
% starting bin, binSize, responseWindow, # bootstrap runs, # baseline window, num k-means runs, number of chance runs.
argList = {... 
    {0.50, 0.10, [2.1 3.1], 0, [0.9 1.9], 10, 1000},...  %The default arguments
    {0.51, 0.10, [2.1 3.1], 0, [0.9 1.9], 10, 1000},...  % Perturb t0
    {0.49, 0.10, [2.1 3.1], 0, [0.9 1.9], 10, 1000},... 
    {0.50, 0.11, [2.1 3.1], 0, [0.9 1.9], 10, 1000},...  % Perturb binSize
    {0.50, 0.09, [2.1 3.1], 0, [0.9 1.9], 10, 1000},... 
          };

M = {};
disp('Computation started...');
tic;
for i = 1:1 %numel(argList) 
  M{i} = {ComputeMetricsForBinaryMixtureTrajectoryClustering(argList{i}{:}), argList{i}};
end
fprintf('Computation completed in %1.3f seconds.\n\n', toc); % 673 seconds.

Mfig       = M{1};
figDir     = GetDataDirForFigure(3);
thisDir    = GetCurrentDirFromPathString(fileparts(mfilename('fullpath')));
outputFile = fullfile(figDir, thisDir, 'recomputedData', 'Mfig.mat');
save(outputFile, 'Mfig');

fprintf('Wrote "%s".\n', outputFile);


