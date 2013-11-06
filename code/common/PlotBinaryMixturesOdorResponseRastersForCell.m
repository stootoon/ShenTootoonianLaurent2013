function PlotBinaryMixturesOdorResponseRastersForCell(whichCell, pnSpt, startTime, endTime, binSize, showRecs, figureId, spikeWidth)
% PlotBinaryMixturesOdorResponseRastersForCell(whichCell, pnSpt, startTime, endTime, binSize, showRecs, figureId, spikeWidth)
%
% Plots a series of rasters showing the response of the selected cell
% to the binary mixtures presented in the binary mixtures experiments,
% as in Figures 2A-D. For a usage example see SpikeTimesBrowser or
% MakeBinaryMixtureRasterExample.
%
% whichCell:  The cell whose responses are to be plotted.
% pnSpt:      The toc matrix of PN spike times, as returned by LoadTocSpikeTimes
% startTime:  The start of the viewing window. Odor onset is at 2 seconds.
% endTime:    The end of the viewing window.   Odor onset is at 2 seconds.
% binSize:    The bin size to use for the reconstructions, in seconds.
% showRecs:   True/false, whether to show the reconstructions.
% figureId:   If not empty, a handle to the figure to plot in.
% spikeWidth: The width of the plotted spikes.
%
% See also: SPIKETIMESBROWSER, MAKEBINARYMIXTURERASTEREXAMPLE.

numTrials = 10;
numOdors  = 27;
numCells  = 168;

pnSptArr = reshape(pnSpt, [], numTrials, numOdors, numCells);

% The mixtures for which reconstructions will be shown.
whichReconcs = [140 30;
                140 80;
                140 120;
                140 140;
                80  140];
if (showRecs)
  % Grab the binned responses for each trial, by passing the list of trials as a cell array.
  whichTrials = arrayfun(@(i) i, 1:numTrials,'UniformOutput',false);
  spkCntsBot = squeeze(CountSpikesInBinsAndAverageAcrossTrials(pnSpt, whichTrials, 1:numOdors, whichCell, 'startTime', startTime, 'endTime', endTime, 'binSize', binSize,'numAllTrials', numTrials,'numAllOdors',numOdors));
end

Q = ComputeSubplotPositions(15,3,[],0.05,0.01,0.1,0.05,0.05,0.01);

whichConc = [140 80]; % Oct Cit
concInd   = GetIndexForBinaryMixtureConcentrationPair(whichConc(1), whichConc(2));

if (isempty(figureId))
  figureId = figure(FindFigureCreate(sprintf('Binary Mixtures: PN %d', whichCell)));
  if (isempty(figureId))
    error('Could not create new figure.');
  end
else
  figure(figureId);
end


clf(figureId);
set(figureId, 'Color',[1 1 1],'Resize','off','NumberTitle','off');

% MATLAB sometimes crashes if it gets interrupted while plotting a
% previous cell's raster, so make sure the figure is not interruptible.
set(figureId, 'Interruptible','off','BusyAction','cancel'); 

ResizeFigure(figureId, 14*10/13,10,'inches');
set(figureId, 'name', sprintf('Binary Mixtures: PN %d (plotting)...', whichCell));

% The citral concentrations to be used in each raster,
% arranged according to their layout in the plot. NAN
% indicates no plot at that location.
citConcs = [nan 000 nan;
            030 030 nan;
            060 060 nan;
            080 080 nan;
            100 100 nan;
            120 120 nan;
            140 140 nan;
            nan nan nan;
            000 140 nan;
            000 140 nan;
            000 140 030;
            000 140 060;
            000 140 080;
            000 140 100;
            nan 140 140;];
citConcs =  citConcs(:);

octConcs = [nan 140 nan;
            000 140 nan;
            000 140 nan;
            000 140 nan;
            000 140 nan;
            000 140 nan;
            000 140 nan;
            nan nan nan;
            140 140 nan;
            120 120 nan;
            100 100 030;
            080 080 060;
            060 060 080;
            030 030 100;
            nan 000 140;];
octConcs =  octConcs(:);

for i = 1:numel(citConcs)
  row = mod(i-1,15)+1;
  col = floor((i-1)/15)+1;
  subplotInd = (row-1)*3+col;

  if (~isnan(citConcs(i)))
    concInd = GetIndexForBinaryMixtureConcentrationPair(octConcs(i), citConcs(i));
    plotSingleRasterInAxis(subplotp(Q,subplotInd), pnSptArr(:,:,concInd,whichCell), startTime, endTime, [octConcs(i) citConcs(i)], spikeWidth);
    % Add concentration labels as necessary
    switch(col)
     case 1
      if (citConcs(i) > 0)
        ylabel(citConcs(i),'Color','g');
      else
        ylabel(octConcs(i),'Color','r');
      end
     case 2
      ylabel({sprintf('oct%d:', octConcs(i)), sprintf('cit%d',citConcs(i))},'FontSize', 9);
     case 3
      ylabel({sprintf('oct%d:', octConcs(i)), sprintf('cit%d',citConcs(i))},'FontSize', 9);
    end
  end

  % Add titles as necessary
  switch col
    case 1
     switch row
      case 1  % First row, insert the panel title
       GhostAxis('axis',subplotp(Q, subplotInd));
       title('A   Component responses','FontSize',12);
      case 2  % First citral row, insert a title
       title('Citral','FontSize',11,'Color','g');
      case 9  % First octanol row, insert a title
       title('Octanol','FontSize',11,'Color','r');
      case 14 % last data row in this column, insert the time axis
       xticks = [-0.5:0.5:2.5];
       xtickLabels = arrayfun(@num2str, xticks,'uniformOutput', false);
       set(gca,'xtick',[-0.5:0.5:2.5]+2,'xticklabel',xtickLabels,'tickLength', [0 0]);
       xlabel('Time(s)','FontSize',10);
     end
   case 2
    if (row == 1)
      title('B   Mixture morph: Oct to Cit','FontSize',12);
    end
   case 3
    switch row
     case 1
      GhostAxis('axis', subplotp(Q,subplotInd));
      title('C   Reconstructions', 'FontSize', 12);
     case 11
      title('D   1:1 Mixture', 'FontSize', 12);
    end
  end

end

% Plot the reconstructions at in the top right.
if (showRecs)
  for i = 1:size(whichReconcs,1)
    row = i+1;
    col = 3;
    subplotInd = (row-1)*3 + col;
    ax = subplotp(Q,subplotInd);
    plotReconstructionInAxis(subplotp(Q,subplotInd), spkCntsBot, startTime, endTime, whichReconcs(i,:), binSize);      
    octConc = whichReconcs(i,1);
    citConc = whichReconcs(i,2);
    ylabel({sprintf('oct%d:', octConc), sprintf('cit%d',citConc)},'FontSize', 9);
  end
end

set(figureId, 'name', sprintf('Binary Mixtures: PN %d', whichCell));

function plotSingleRasterInAxis(ax,spt,t0,t1,octCitMix, spikeWidth)
if (any(isnan(octCitMix)))
  return;
end

octConc = octCitMix(1);
citConc = octCitMix(2);

set(gcf,'CurrentAxes',ax);
patchColor = [0.902 0.902 0.9804]; % lavender

citralColor30 = [0.85 1.0 0.85];
citralColor140= [0.0 1.0 0.0];

octanalColor30 = [1.0 0.85 0.85];
octanalColor140= [1.0 0.0 0.0];

octanalColor = interp1([0;30;140], [1 1 1; octanalColor30; octanalColor140], octCitMix(1));
citralColor  = interp1([0;30;140], [1 1 1;  citralColor30;  citralColor140], octCitMix(2));

cla;
patch([2 2.3 2.3 2],    [0.25 0.25 10 10],   patchColor, 'EdgeColor', 'none'); % The odor pach

% Spt has the spike times in its rows and the columns are trials. 
% Add the columns as a complex number to make life easier for us.

spt    = bsxfun(@plus, spt, (0:10-1)*sqrt(-1));
spt    = spt(:);
indVld = find(real(spt)>=t0+0.3 & real(spt)<=t1); % +0.3 to allow some space for the mixture indicator

spikeTime = real(spt(indVld));
spikeTrial= imag(spt(indVld));


spikeHeight= 1;

% Plot with lines instead of patches, looks better
% spikeWidth = 0.03;
% X = bsxfun(@plus, spikeTime, [0 1 1 0]*spikeWidth);
% Y = bsxfun(@plus, spikeTrial, [0 0 1 1]*spikeHeight);
% patch(X',Y',[0 0 0],'EdgeColor','none');

X = spikeTime(:)*[1 1 nan]; 
Y = bsxfun(@plus, spikeTrial(:), [0 spikeHeight nan]);
X = X';
Y = Y';
hold on;
plot(X(:),Y(:),'k','LineWidth', spikeWidth);
patch([0.025 0.3 0.3 0.025]+t0, 5-0.95*[0 0 5 5]*citConc/140,    citralColor, 'EdgeColor', 'none');
patch([0.025 0.3 0.3 0.025]+t0, 0.95*[0 0 5 5]*octConc/140+5,   octanalColor, 'EdgeColor', 'none');

% Setup the xticks to span [t0 - t1]
xt0 = when(t0<0, -ceil(abs(t0)/0.5), floor(t0/0.5));
xt1 = when(t1<0, -ceil(abs(t1)/0.5), floor(t1/0.5));
xtick = [xt0:xt1]*0.5;
set(gca,'xtick', xtick, 'xticklabel',[],'ytick',[],'yticklabel',[],'ticklength',[0 0]);

grid on;
box on;
xlim([t0,t1]);
axis ij;

function plotReconstructionInAxis(ax, cntsBot, t0, t1, octCitMix, binSize)
% Get the indices for the components and the mixture
octConc = octCitMix(1);
citConc = octCitMix(2);

citralColor30 = [0.85 1.0 0.85];
citralColor140= [0.0 1.0 0.0];

octanalColor30 = [1.0 0.85 0.85];
octanalColor140= [1.0 0.0 0.0];

octanalColor = interp1([0;30;140], [1 1 1; octanalColor30; octanalColor140], octCitMix(1));
citralColor  = interp1([0;30;140], [1 1 1;  citralColor30;  citralColor140], octCitMix(2));

citInd  = GetIndexForBinaryMixtureConcentrationPair(0, citConc);
octInd  = GetIndexForBinaryMixtureConcentrationPair(octConc, 0);
mixInd  = GetIndexForBinaryMixtureConcentrationPair(octConc, citConc);

mixCnt  = mean(cntsBot(:,mixInd,:),3);
recCnt  = mean(sum(cntsBot(:,[citInd octInd],:),2),3);

numBins = size(cntsBot,1);
t = (0:numBins-1)*binSize + t0;
set(gcf,'CurrentAxes',ax);

cla;
patchColor = name2rgb('lavender');
ymax = 2.1;
patch([2 2.3 2.3 2],    0.0+0.975*[0 0 1 1]*ymax, patchColor, 'EdgeColor', 'none'); hold on;
patch([0.025 0.3 0.3 0.025]+t0, ymax/2 + 0.95*[0 0 ymax ymax]/2*citConc/140,   citralColor, 'EdgeColor', 'none');
patch([0.025 0.3 0.3 0.025]+t0, ymax/2 - 0.95*[0 0 1 1]*ymax/2*octConc/140,   octanalColor, 'EdgeColor', 'none');

%plot(t(t>t0+0.3), mixCnt(t>t0+0.3), 'Color','b');  hold on;
set(area(t(t>t0+0.3), mixCnt(t>t0+0.3)),'EdgeColor','none','FaceColor', name2rgb('gray50')); 
hold on;
plot(t(t>t0+0.3), recCnt(t>t0+0.3), 'Color','k');

% Setup the xticks to span the xlim
xt0 = when(t0<0, -ceil(abs(t0)/0.5), floor(t0/0.5));
xt1 = when(t1<0, -ceil(abs(t1)/0.5), floor(t1/0.5));
xtick = [xt0:xt1]*0.5;
set(gca,'xtick',xtick,'xticklabel',[],'ytick',[],'yticklabel',[],'ticklength',[0 0]);
grid on;
box on;
ylim([0 2.1]);
xlim([t0 t1]);
