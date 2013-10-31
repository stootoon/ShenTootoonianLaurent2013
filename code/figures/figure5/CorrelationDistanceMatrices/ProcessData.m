function ProcessData()
% function ProcessData()

figDir    = GetDataDirForFigure(5);
currDir   = GetCurrentDirFromPathString(fileparts(mfilename('fullpath')));
targetDir = fullfile(figDir, currDir, 'recomputedData');

fileName   = fullfile(targetDir, 'corrDistResultsSingleComponents.mat');
startTime  = tic;
whichOdors = {'A','B','C','D','W','X','Y','Z'};
computeCorrelationDistancesForOdors(whichOdors, fileName);
fprintf('Wrote "%s" in %1.3f seconds.\n', fileName, toc(startTime));

fileName   = fullfile(targetDir, 'corrDistResultsWseries.mat');
startTime  = tic;
whichOdors = {'W','WX','WXY','WXYZ','AWXYZ'};
computeCorrelationDistancesForOdors(whichOdors, fileName);
fprintf('Wrote "%s" in %1.3f seconds.\n', fileName, toc(startTime));

fileName   = fullfile(targetDir, 'corrDistResultsWfamilyAfamily.mat');
startTime  = tic;
whichOdors = {'W','X','Y','Z','WX','WY','WZ','XY','XZ','YZ','WXY','WYZ','WXYZ',...
              'A','B','C','D','AB','AC','AD','BC','ABC','ACD','ABCD'};
computeCorrelationDistancesForOdors(whichOdors, fileName);
fprintf('Wrote "%s" in %1.3f seconds.\n', fileName, toc(startTime));

fileName   = fullfile(targetDir, 'corrDistResultsPartialOverlaps.mat');
startTime  = tic;
whichOdors = {'WXYZ','DWYZ','ABCD','ABCX','BCWX','BDWX'};
computeCorrelationDistancesForOdors(whichOdors, fileName);
fprintf('Wrote "%s" in %1.3f seconds.\n', fileName, toc(startTime));

fileName   = fullfile(targetDir, 'corrDistResultsTransientOverlaps.mat');
startTime  = tic;
whichOdors = {'D','AD','ACD','ABCD','Z','WZ','WYZ','DWYZ'};
computeCorrelationDistancesForOdors(whichOdors, fileName);
fprintf('Wrote "%s" in %1.3f seconds.\n', fileName, toc(startTime));

function computeCorrelationDistancesForOdors(whichOdors, fileName)
whichOdorInds = cellfun(@(x) odor_name_to_index(['odor' x]), whichOdors)
X  = LoadTocSpikeTimes('rawpn');
Y1 = CountSpikesInBinsAndAverageAcrossTrials(X, 2:4, whichOdorInds, 1:174, 'startTime', 2, 'endTime', 5);
Y2 = CountSpikesInBinsAndAverageAcrossTrials(X, 5:7, whichOdorInds, 1:174, 'startTime', 2, 'endTime', 5);

f  = @(M) min(M(:));
Y1 = reshape(Y1,174,[],numel(whichOdorInds));
Y2 = reshape(Y2,174,[],numel(whichOdorInds));
D  = ComputeCorrelationDistance(Y1,Y2);
Ds = SummarizeCorrelationDistanceMatrices(D,f);

save(fileName, 'D', 'Ds', 'whichOdors');