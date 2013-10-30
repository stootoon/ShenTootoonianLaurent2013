function [r2, e2, snr, signalPower, noisePower] = ComputeLinearityIndicesFromBayes1LaplaceResults2(Results, pnCbotDataFile, varargin)
% [r2, e2, snr, signalPower, noisePower] = ComputeLinearityIndicesFromBayes1LaplaceResults(Results, pnCbotDataFile, varargin)
% 
% Compute statistics measuring the quality of the fits of the mixture
% responses to the componentes.  The statistics computed are the R^2
% (SSR/SST) and E^2 (SSE/SST) of the fits and the SNRs of the
% responses. The R^2 and E^2 values are read out direclty from the
% results file, the SNR is computed here. The results are returned in
% the NUMCELLS x NUMMIXTURES matrices R2, E2, and SNR. 
%
% SIGNALPOWER and NOISEPOWER are (NUMCELLS x NUMIXTURES) x 3 matrices
% containing the signal power and the noise power computed for each of
% the components and the mixture, in a particular response. For
% example if the i'th row contains the data for the response of cell j
% to mixture k, signalPower(i,1) contains the signal power for that
% cell in response to the octanol component of that mixture,
% signalPower(i,2) the response to the citral component, and
% signalPower(i,3) to the mixture. Similarly for noisePower.

p = inputParser;
p.addOptional('tbl',   [0 2; 8 13]); % Baseline time bins
p.addOptional('tresp', [2 4]); % Baseline time bins 

% Summary functions will be applied to a three-element vector
% containing data for each of the component odors and the mixture
% odor. The ratio of the signal summary to the noise summary yields the SNR.
p.addOptional('noiseSummaryFunction', @mean);        % noiseSummary:  average noise across all three conditions.
p.addOptional('signalSummaryFunction', @(x) x(end)); % signalSummary: the signal during the mixture condition.
p.KeepUnmatched = true;
p.parse(varargin{:});

noiseSummaryFunction = p.Results.noiseSummaryFunction;
signalSummaryFunction= p.Results.signalSummaryFunction;

tbl   = p.Results.tbl;
tresp = p.Results.tresp;

fitOptions = UnpackStructureFieldsAsNameValuePairs(p.Unmatched);

if (size(Results,2)~=1)
  error('Expected results to be a column cell array.');
end

numResults = numel(Results);

r2 = zeros(size(Results));
e2 = zeros(size(Results));
snr= zeros(size(Results));

[whichMixtures, whichCells] = meshgrid(1:16, 1:168); % 16 is the total number of binary mixture series.
[mixtureInds, mixtureVals]  = GetBinaryMixtureIndsForResponseLinearity;

whichMixtures = whichMixtures(:);
whichCells    = whichCells(:);

Data = load(pnCbotDataFile);

% Figure out the baseline and response bins
indBl = [];
for i = 1:size(tbl,1)
  indBl = [indBl find(Data.binStarts>=tbl(i,1)  & Data.binStarts<tbl(i,2))];
end
indResp     = find(Data.binStarts>=tresp(1) & Data.binStarts<tresp(2));

% Now read out the R2 and E2 values from the input file, and also
% compute the SNRs.

noisePower  = zeros(numResults, 3); % 3 = 2 components + 1 mixture
signalPower = zeros(numResults, 3); 
for i = 1:numResults
  thisMixture = mixtureVals(whichMixtures(i),:);
  cmpInds = [GetIndexForBinaryMixtureConcentrationPair(thisMixture(1),0);...
             GetIndexForBinaryMixtureConcentrationPair(0,thisMixture(2))]';
  respInd = GetIndexForBinaryMixtureConcentrationPair(thisMixture(1), thisMixture(2));

  % Xybo is a numBins x 3 matrix containing the cell's mean response
  % across trials to the two components and the mixture, in each time
  % bin. 
  Xybo        = squeeze(mean(Data.pnCbot(whichCells(i), :, [cmpInds respInd], :),4));
  muNoise     = mean(Xybo(indBl,:),1);
  noisePower(i, :) =  var(Xybo(indBl,:),[],1);  
  signalPower(i,:) = sum(bsxfun(@minus, Xybo(indResp, :), muNoise).^2)/numel(indResp); % Mean squared deviation from the average baseline response.
  
  noise = noiseSummaryFunction(noisePower(i,:));
  signal= signalSummaryFunction(signalPower(i,:));
    
  snr(i) = signal/noise;
  r2(i)  = Results{i}.r2;
  e2(i)  = Results{i}.e2;
end

snr   = reshape(snr,   168, 16);
r2    = reshape(r2,    168, 16);
e2    = reshape(e2,    168, 16);
