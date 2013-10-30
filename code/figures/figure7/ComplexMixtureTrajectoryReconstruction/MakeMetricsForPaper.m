function MakeMetricsForPaper(varargin)
% function MakeMetricsForPaper(['dataDir' = 'originalData'])
p = inputParser;
p.addOptional('dataDir', 'originalData');
p.parse(varargin{:});

figDir  = GetDataDirForFigure(7);
dataDir = p.Results.dataDir;
currDir = GetCurrentDirFromPathString(fileparts(mfilename('fullpath')));

data           = load(fullfile(figDir, currDir, dataDir, 'dataForPnRec.mat'));
unshuffledSsen = data.ssenData{1};
shuffledSsen   = data.ssenData{3};
fprintf('Median of SSE/SST for unshuffled data: %1.3f.\n', median(unshuffledSsen));
fprintf('Median of SSE/SST for   shuffled data: %1.3f.\n', median(shuffledSsen));
fprintf('p-value for median of SSE/SST for shuffled data being < for unshuffled: %1.3e\n', ranksum(unshuffledSsen, shuffledSsen)/2); %/2 to get sided p-value.