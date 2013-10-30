function [snr, SignalPower, NoisePower] = ComputeSnrForAllCellsAndMixtures(varargin)
% [snr, SignalPower, NoisePower] = ComputeSnrForAllCellsAndMixtures(varargin)
%
% Computes the SNR for each mixture response of each single PN. This
% is done by examining, for each PN and each mixture response, the
% single component responses and the mixture response. The mean and
% variance of the power in the baseline is then computed, to yield a
% 'noisePower' signal, containing one value for each component and the
% mixture. We then compute the variance of the the PN during the
% response window, but using the baseline mean - hence the mean of the
% squared deviation from the baseline mean, to yield the 'signal
% power' vecotr. A summary function is then applied to the vector or
% noise responses, and a second summary function to the signal power
% vector, each yielding a single number: the signal and the noise. We
% then divide the signal by the noise to get the SNR. 
%
% INPUTS: (can be overriden via name-value pairs):
%
%  dataFile([pnComplexMixtures0_0to13_0_50ms.mat]): The file containing
%  the PN data.
%   
%  baselineIntervals[{[0 2],[8 13]}]: The (possibly several) time
%  intervals to use for computing the baseline mean and variance.
%
%  responseInterval[ [2,4] ]: The time interval defining the response
%  window.
%   
%  noiseSummaryFunction[@mean]: The function to apply to the vector of
%  noise powers computed for each component and the mixture to yield a
%  single, 'noise' value.
%
%  signalSummaryFunction[@,ax]: The function to apply to the vector of
%  signal powers computed for each component and the mixture to yield a
%  single, 'signal' value.
%
% OUTPUTS:
% 
%  SNR: A NUMPNS x NUMMIXTURES matrix of computed SNR values
%
%  SignalPower: A NUMPNS x NUMMIXTURES cell matrix of signal vectors
%  computed for the component and mixture responses.
%
%  NoisePower:  A NUMPNS x NUMMIXTURES cell matrix of noise vectors
%  computed for the component and mixture responses.

p = inputParser;
p.addOptional('dataFile','pnComplexMixtures0_0to13_0_50ms.mat');
p.addOptional('baselineIntervals',{[0 2],[8 13]});
p.addOptional('responseInterval',[2 4]);
p.addOptional('noiseSummaryFunction',  @mean);
p.addOptional('signalSummaryFunction', @(x) x(end));
p.parse(varargin{:});

Data = load(p.Results.dataFile);
baselineIntervals = p.Results.baselineIntervals;
responseInterval  = p.Results.responseInterval;
noiseSummaryFunction  = p.Results.noiseSummaryFunction;
signalSummaryFunction = p.Results.signalSummaryFunction;

responseBins = find(Data.tall>=responseInterval(1) & Data.tall<responseInterval(2));
baselineBins = cellfun(@(x) Columnize(find(Data.tall>=x(1) & Data.tall<x(2))), baselineIntervals, 'UniformOutput', false);
baselineBins = vertcat(baselineBins{:});

cmpsInOdors = GetOdorNamesAsBinaryVectors('full');

pnCbot = Data.pnCbot;
[numCells, numBins, numOdors, numTrials] = size(pnCbot);
Xm = squeeze(mean(pnCbot,4));

mixtureOdors = 13:44; % AB to ABCDWXYZ
numMixtureOdors = mixtureOdors(end)-mixtureOdors(1)+1;

SignalPower = cell(numCells, numMixtureOdors); 
NoisePower  = SignalPower;

snr = zeros(numCells, numMixtureOdors);
for i = 1:numCells
  for j = mixtureOdors % Mixtures
    cmpsInThisOdor = find(cmpsInOdors(j,:));
    indCmps = cmpsInThisOdor+3;
    Xybo = squeeze(Xm(i,:,[indCmps j]));
    jj = j - mixtureOdors(1) + 1;
    [NoisePower{i,jj}, muNoise] = computeNoisePower(Xybo);
    SignalPower{i,jj} = computeSignalPower(Xybo, muNoise);    
  end
end

np = cellfun(noiseSummaryFunction,  NoisePower);
sp = cellfun(signalSummaryFunction, SignalPower);
snr = sp./np;

function signalPower = computeSignalPower(Xybo, muNoise)
XyResp      = Xybo(responseBins,:);
signalPower = sum(bsxfun(@minus, XyResp, muNoise).^2)/numel(responseBins);
end

function [noisePower, muNoise] = computeNoisePower(Xybo)
XyNoise    = Xybo(baselineBins,:);
muNoise    = mean(XyNoise,1); 
noisePower = var(XyNoise);
end

end


