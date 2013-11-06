function MakeBinaryMixtureRasterExample()
% function MakeBinaryMixtureRasterExample()
%
% This function makes the figure showing the raster for a binary
% mixture response (Figure 2A-D). 

whichCell  = 128; % The cell whose data will be plotted.
startTime  = 1.2;   % Odor onset is at t = 2.0
endTime    = 4.5;
binSize    = 0.05;  % For the reconstructions
showRecs   = true;
spikeWidth = 1;

pnSpt      = ConvertSpikeTimesFromSparseToFull(LoadTocSpikeTimes('rawpn_binary_mixtures'));
PlotBinaryMixturesOdorResponseRastersForCell(whichCell, pnSpt, startTime, endTime, binSize, showRecs, [], spikeWidth);
set(gcf, 'Name', 'Figure 2A-D: Example PN Binary Mixture Response Rasters'); 
