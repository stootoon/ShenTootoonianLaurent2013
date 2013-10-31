function ProcessData()
% function ProcessData()

figDir    = GetDataDirForFigure(7);
currDir   = GetCurrentDirFromPathString(fileparts(mfilename('fullpath')));
targetDir = fullfile(figDir, currDir, 'recomputedData');

startTime = 2;
endTime   = 5;
ml        = OdorGroupingFunction1(GetOdorsList);
respTh    = 4; % Number of responsive tirals for a cell to be deemed responsive

startTime = tic;

pnSpt     = ConvertSpikeTimesFromSparseToFull(LoadTocSpikeTimes('rawpn'));
pnCbot50  = CountSpikesInBinsAndAverageAcrossTrials(pnSpt,{1,2,3,4,5,6,7},1:44,1:174,'startTime',2,'endTime',5,'binSize',0.05);
pnCbot100 = CountSpikesInBinsAndAverageAcrossTrials(pnSpt,{1,2,3,4,5,6,7},1:44,1:174,'startTime',2,'endTime',5,'binSize',0.10);

kcSpt     = ConvertSpikeTimesFromSparseToFull(LoadTocSpikeTimes('rawkc'));
kcCbot50  = CountSpikesInBinsAndAverageAcrossTrials(kcSpt,{1,2,3,4,5,6,7},1:44,1:209,'startTime',2,'endTime',5,'binSize',0.05);
kcCbot100 = CountSpikesInBinsAndAverageAcrossTrials(kcSpt,{1,2,3,4,5,6,7},1:44,1:209,'startTime',2,'endTime',5,'binSize',0.10);

% Worker functions
MFRVSML    = @(cbot, whichMl) squeeze(mean(mean(cbot(:,:,ml==whichMl,:),   4)));          % First average over trials, then over cells
RESPONSIVE = @(cbot, whichMl) squeeze(mean( sum(cbot(:,:,ml==whichMl,:)>0, 4)>=respTh)); % Compute responsivity using the trials, then average over cells.
SILENT     = @(cbot, whichMl) squeeze(mean( sum(cbot(:,:,ml==whichMl,:)>0, 4)==0));      % Compute silence using the trials, then average over cells.
CUMRESP    = @(cbot, whichMl) mean(Columnize(sum(sum(cbot(:,:,ml==whichMl,:)>0,4)>=respTh, 2)>0)); % Averages the overall responsivity level at a given ML

%% Build the structure
Data = struct;

Data.mixtureLevels = unique(ml);
Data.respTh        = respTh;

Data.startTime = startTime;
Data.endTime   = endTime;

% Actually do the computation.
Data.pnSz50       = size(pnCbot50);
Data.pnSz100      = size(pnCbot100);
Data.pnMfrvsml    = arrayfun(@(i) MFRVSML(pnCbot50,    i), unique(ml), 'UniformOutput', false);
Data.pnResponsive = arrayfun(@(i) RESPONSIVE(pnCbot50, i), unique(ml), 'UniformOutput', false);
Data.pnSilent     = arrayfun(@(i) SILENT(pnCbot100,    i), unique(ml), 'UniformOutput', false);
Data.pnCumResp    = arrayfun(@(i) CUMRESP(pnCbot50,    i), unique(ml));

Data.kcSz50       = size(kcCbot50);
Data.kcSz100      = size(kcCbot100);
Data.kcMfrvsml    = arrayfun(@(i) MFRVSML(kcCbot50,    i), unique(ml), 'UniformOutput', false);
Data.kcResponsive = arrayfun(@(i) RESPONSIVE(kcCbot50, i), unique(ml), 'UniformOutput', false);
Data.kcSilent     = arrayfun(@(i) SILENT(kcCbot100,    i), unique(ml), 'UniformOutput', false);
Data.kcCumResp    = arrayfun(@(i) CUMRESP(kcCbot50,    i), unique(ml));

Data.mfrvsmlFun    = MFRVSML;
Data.responsiveFun = RESPONSIVE;
Data.silentFun     = SILENT;
Data.dateComputed  = datestr(now);

Data.notes = 'Data for time course of PN and KC responsivity.\nA cell is considered RESPONSIVE if it produces APs in at least RESPTH of the 7 trials avaialable in each time bin.\nIt''s considered SILENT if it produces 0 APs in that bin.\nBin size is 50 ms, but 100 ms for computing the fraction of SILENT cells.\nThe fields of DATA contain the results computed for PNs and KCs.\nEach cell corresponds to a mixture level, ordered according to the MIXTURELEVELS field.\nEach of the elements of each cell array contains the results for each of the odors in that ML.\n\n';

outputFile = fullfile(targetDir, 'dataPnKcResponseTimeCourse.mat');
save(outputFile, 'Data');
fprintf('Wrote "%s" in %1.3f seconds.\n', outputFile, toc(startTime));