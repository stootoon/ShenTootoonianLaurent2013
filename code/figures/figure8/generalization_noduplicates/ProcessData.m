function ProcessData()
% function ProcessData()

global targetDir
figDir    = GetDataDirForFigure(8);
currDir   = GetCurrentDirFromPathString(fileparts(mfilename('fullpath')));
targetDir = fullfile(figDir, currDir, 'recomputedData');

numMpCores = 6;
if (MatlabPoolWrapper('size')>0)
  MatlabPoolWrapper('close');
end
MatlabPoolWrapper('open', numMpCores);

disp('Preparing data...');
startTime = tic;
dataFile = PrepareDatasets(); % ~24.7 seconds
fprintf('Wrote "%s" in %1.1f secs.\n', dataFile, toc(startTime));

disp('Computer per-component categorization performance...');
startTime = tic;
catResultsFile = ParComputePerComponentCategorizationResults(dataFile); % ~23000 secs
fprintf('Wrote "%s" in %1.1f secs.\n', catResultsFile, toc(startTime));

disp('Computer per-component generalization performance...');
startTime = tic;
genResultsFile = ParComputePerComponentGeneralizationResults(dataFile); % ~3200 secs
fprintf('Wrote "%s" in %1.1f secs.\n', genResultsFile, toc(startTime));

function outputFile = PrepareDatasets()
global targetDir

kcSpt = ConvertSpikeTimesFromSparseToFull(LoadTocSpikeTimes('rawkc'));
pnSpt = ConvertSpikeTimesFromSparseToFull(LoadTocSpikeTimes('rawpn'));

disp('Creating binned KC spike times...');
[kcraw,kcsw] = createDataset(kcSpt);

disp('Creating binned PN spike times...');
[pnraw,pnsw] = createDataset(pnSpt);

notes = 'pnraw and kcraw are the raw pn and kc spike times, binned in 25 msec windows between 2.0 and 4.1. The dimensions are cells, bins, odors, trials. pnsw and kcsw are sliding windows on this data, where the spikes in a 100 msec bin have been computed as the bin slides forward in 25 msec steps. The dimension meanings are the same, but there are 3 fewer bins due the use of a 100 msec window.';

disp('Saving to binnedSpikeTimes.mat...');

outputFile = fullfile(targetDir, 'binnedSpikeTimes.mat');
save(outputFile, 'kcraw', 'kcsw', 'pnraw', 'pnsw', 'notes');

function [raw, sw] = createDataset(spt)
numCells = size(spt,2)/308;
raw = CountSpikesInBinsAndAverageAcrossTrials(spt,map(@Identity,1:7),1:44,1:numCells,'startTime',2.0,'endTime',4.1,'binSize',0.025);
sw = 0*raw;
sw = sw(:,1:end-3,:,:);
for i = 1:size(sw,2)
  sw(:,i,:,:) = sum(raw(:,i:i+3,:,:),2);
end

function outputFile = ParComputePerComponentCategorizationResults(dataFile)
global targetDir

numReps    = 50;
lambda     = 1; 
verbose    = 0;
CatResults = {};
cmps       = 'ABCDWXYZ';
whichCmps  = 1:numel(cmps);

load(dataFile);

for i = 1:numel(whichCmps)
  icmp = whichCmps(i);
  SetRandomSeed1(icmp);
  cmp  = cmps(icmp);
  Pn   = ComputeCategorizationPerformanceForComponent(pnsw, cmp, numReps, lambda, verbose);
  Kc   = ComputeCategorizationPerformanceForComponent(kcsw, cmp, numReps, lambda, verbose);
  CatResults{i} = {cmp, Pn, Kc};  
end
CatResults = CatResults(:);

outputFile = fullfile(targetDir, 'CatResults.mat');
save(outputFile, 'CatResults');

function outputFile = ParComputePerComponentGeneralizationResults(dataFile)
global targetDir

numReps    = 50;
lambda     = 1; 
verbose    = 0;
GenResults = {};
cmps       = 'ABCDWXYZ';
whichCmps  = 1:numel(cmps);

load(dataFile);

for i = 1:numel(whichCmps)
  icmp = whichCmps(i);
  SetRandomSeed1(icmp);
  cmp  = cmps(icmp);
  Pn   = ComputeGeneralizationPerformanceForComponent(pnsw, cmp, numReps, lambda, verbose);
  Kc   = ComputeGeneralizationPerformanceForComponent(kcsw, cmp, numReps, lambda, verbose);
  GenResults{i} = {cmp, Pn, Kc};  
end
GenResults = GenResults(:);

outputFile = fullfile(targetDir, 'GenResults.mat');
save(outputFile, 'GenResults');
