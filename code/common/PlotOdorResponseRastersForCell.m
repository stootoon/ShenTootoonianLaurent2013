function ax = PlotOdorResponseRastersForCell(tocSpikeTimes, whichCell, timeWindow, odorToHighlight, varargin)
% ax = PlotOdorResponseRastersForCell(spt, whichCell, timeWindow, odorToHighlight, ...)
%
% Plots the spike times of a single cell as an array of rasters,
% arranged Kai's way. Adapted from Kai's code.
%
% tocSpikeTimes: Dataset to uset.
% 
% whichCell: Which cell to plot
%
% timeWindow: The window that the plots should cover, e.g. [0 1].
%
% odorToHighlight: The colors of all the rasters for odors
% containing this component will be tinted red. 
%
% 'ax' is a vector of handles to the subplots containing the rasters.

opts = inputParser;
opts.addOptional('odorOnsetTime', 2);
opts.addOptional('spikeColor',  [0 0 0]);
opts.addOptional('spikeWidth',  1);
opts.addOptional('clearFigure', true);
opts.addOptional('figureId', []);
opts.parse(varargin{:});

opts = opts.Results;

odorStrs = GetOdorsList();
numOdors = numel(odorStrs);
numTrialsPerCond = 7;  

S             = tocSpikeTimes(:,(1:308)+(whichCell-1)*308);
odorOnsetTime = opts.odorOnsetTime;
stimDur       = 0.5;

%% Set position of different odors
iRows = [1 3 7 2 5 4 6 8 10 9 11 12   3 4 2 1 5 6 12 8 7 10 11 9    3 2 1 5 6 10 8 9        2 3 5 6 9 8         2 3 9 8 5          6];
iCols = [1*ones(1,12)                 2*ones(1,12)                  3*ones(1,8)             4*ones(1,6)         5*ones(1,5)        6];

nCols = max(iCols);
nRows = max(iRows);

vSpacing    = 0.03;
hSpacing    = 0.005;
hMargin     = 0.04;
vMargin     = 0.04;

plotWidth   = (1-2*hMargin - (nCols-1)*hSpacing)/nCols;
plotHeight  = (1-2*vMargin - (nRows-1)*vSpacing)/nRows; 

plotLeft   = (iCols-1)*plotWidth+(iCols>1).*(iCols-1)*hSpacing + hMargin;
plotBottom = 1 - vMargin - (iRows*plotHeight + (iRows>1).*(iRows-1)*vSpacing);

if (~isempty(opts.figureId))
  sfigure(opts.figureId);
end

if (opts.clearFigure)
  clf;
end
set(gcf,'Interruptible','off','BusyAction','cancel');

ff = gcf;
trialInd = 1;

% Convert S to a cell array if necessary
C = cell(1,numTrialsPerCond*numOdors);
for i = 1:size(S,2)
  ind = find(S(:,i));
  C{i} = S(ind,i)';
end
S = C;

ax = zeros(1,length(odorStrs));
for oCtr = 1:length(odorStrs)
  oStr = odorStrs{oCtr};
  oStr = strrep(oStr, 'high' ,'  4 x');
  oStr = strrep(oStr, 'ALL', 'ABCDWXYZ');
  ax(oCtr) = subplot('Position', [plotLeft(oCtr), plotBottom(oCtr), plotWidth, plotHeight]);

  % Overlap/Gap patch
  patchY = [0.2 numTrialsPerCond numTrialsPerCond 0.2];
  patchX = odorOnsetTime+[0 0 stimDur stimDur]; 
  transp = 0.2;

  if (opts.clearFigure)
    p1 = patch(patchX, patchY, [0.8 0.8 0.8]); 
    set(p1, 'EdgeColor', 'none'); 
  end

  hold on;
  X = [];
  Y = [];
  for j = 1:7
    if (~isempty(S{trialInd}))
      spikeTimes = S{trialInd};
      X = [X spikeTimes];
      Y = [Y (j-1)*ones(size(spikeTimes))];
    end
    trialInd = trialInd+1;
  end
  X = X(:)';
  X = [X;X;nan*X];
  X = X(:)';
  Y = Y(:)';
  Y = [Y;Y+1;nan*Y];
  Y = Y(:)';
  pp = plot(X,Y, 'Color', opts.spikeColor, 'LineWidth', opts.spikeWidth);

  axis([timeWindow(1) timeWindow(2) 0 numTrialsPerCond]);
  axis ij;
    
  % Add text, patches of color etc.
  hText = text(timeWindow(1), 0, oStr(5:length(oStr)),'VerticalAlignment', 'bottom'); 
  set(hText,'FontName','Arial');
  set(hText, 'FontSize', 12);
  
  % Set color according to odor component
  if (oCtr == 44)
    newStr = 'ABCDWXYZ';
  else
    newStr = oStr(5:length(oStr));
  end
  
  indCmpInOdor = findstr(newStr,odorToHighlight);
  
  if isempty(indCmpInOdor)
    set(hText, 'Color', 'k');
    set(hText, 'FontWeight', 'bold');
  else
    set(hText, 'Color', name2rgb('magenta'));
    set(hText, 'FontWeight', 'bold');
  end
    
  h = gca;
  set(h, 'XTick', -2:0.5:10);
  if (oCtr~=44)
    set(h, 'XTickLabel', []);
  else
    set(h, 'XTickLabel', arrayfun(@(t) num2str(t - odorOnsetTime), get(gca,'xtick'), 'UniformOutput', false));
  end
  set(h, 'XMinorTick', 'off');

  set(h, 'YTick', [], 'XColor', [0 0 0], 'YColor', [0 0 0],'Box','on', 'TickLength', [0 0]);
end
xl = xlabel('Time (s)');
set(xl,'FontSize',12);
set(ff,'color','white');
orient landscape;

