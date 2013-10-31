function PlotLinearityVsSnrAsClickMapForBayes1Laplace(Results, r2, snr, varargin)
% PlotLinearityVsSnrAsClickMapForBayes1Laplace(Results, r2, snr, varargin)
%
% Given the NUMCELLS X NUMMIXTURES matrix of structures RESULTSNL
% (e.g. as returned by COMPUTELINEARITYINDICESFROMRESULTS), plots the
% points as SNR vs R2, as a click map so that clicking a point brings
% up the associated response. 

p = inputParser;
p.addOptional('pvals',[]);
p.addOptional('pnCbotDataFile','delayRegressData3.mat');
p.KeepUnmatched = true;
p.parse(varargin{:});
pvals = p.Results.pvals;

fitOptions = UnpackStructureFieldsAsNameValuePairs(p.Unmatched);
pnCbotDataFile = p.Results.pnCbotDataFile;

[numCells, numMixtures] = size(Results);
bm = cellfun(@(x) x.bestModel, Results);
snrDb = 10*log10(abs(snr));

plotMode = 'standard';
hHighlightedCell = [];
hHighlightedCellMix = [];
clickMap = makePlot();

function clickMap = makePlot()
clickMap = sfigure(FindFigureCreate('R2 vs SNR Clickmap'));
clf; whitebg(gcf,'w'); set(gcf,'Tag','R2SnrClickmap');
rnd = rand(size(r2(:)));
if (isempty(pvals) | ~isequal(plotMode, 'pval'))
  cols = [rnd(:).*(bm(:)==1) 0*rnd(:) rnd(:).*(bm(:)~=1)]*0.5;
  cols(bm==1,:) = bsxfun(@plus, [0.5 0 0], cols(bm==1,:));
  cols(bm>1,:) = bsxfun(@plus, [0 0 0.5], cols(bm>1,:));
else
  cols = interp1([0 0.1 0.2 0.5 1],[0.5 0 0; 1 0 0; 1 1 1; 0 0 1; 0 0 0.5],pvals(:));
end
h = scatter3(snrDb(:), r2(:), randn(size(r2(:))),20,cols,'filled');
set(h,'hittest','off');
set(gca,'ButtonDownFcn', @(h,evt) mapPointToPlot(get(gca,'CurrentPoint')));
set(gcf, 'KeypressFcn', @keypressFcn);
view(2);
axis square;
ylim([-0.05 max(ylim)]);
xlabel('SNR (dB)');
ylabel('R^2');

switch(plotMode)
 case 'standard'
  title('R^2 vs SNR(dB), colored by CELL. Click on a point to see the corresponding response');
 case 'pval'
  title('R^2 vs SNR(dB), colored by R^2 pvalue. Click on a point to see the corresponding response');
 otherwise
  title('Click on a point to see the corresponding response.');
end
TightenAxesToFigure;
end

function keypressFcn(src, evt)
if (isequal(evt.Character,' '))
  switch plotMode
   case 'standard'
    if (isempty(pvals))
      disp('No pvalues provided.');
    else
      disp('Plotting in pval mode.');
      plotMode = 'pval';
      makePlot();
    end
   case 'pval'
    disp('Plotting in standard mode.');
    plotMode = 'standard';
    makePlot();
   otherwise
    error('Unknown plot mode "%s".\n', plotMode);
  end
end
end 

function ind = nearestInd(x,y)
dx = diff(xlim);
dy = diff(ylim);
ind = argmin((snrDb(:)-x).^2/dx^2+(r2(:)-y).^2/dy^2,1);
fprintf('Nearest ind to (%1.3f %1.3f) is %d: (%1.3f, %1.3f).\n', x,y,ind, snrDb(ind), r2(ind));
end

function imix = mixFromInd(i)
imix = floor((i-1)/numCells)+1;
end

function icell = cellFromInd(i)
icell = mod((i-1), numCells)+1;
end

function mapPointToPlot(pt)
nearestInd = nearestInd(pt(1),pt(3));
cellInd    = cellFromInd(nearestInd);
mixInd     = mixFromInd(nearestInd);
fprintf('Clicked (%1.3f, %1.3f): cell %d, mixture %d.\n',pt(1),pt(3),cellInd, mixInd);
ShowModelFitsBayes1Laplace(cellInd, mixInd, Results{cellInd, mixInd}.bestModel,pnCbotDataFile, fitOptions{:});
highlightCell(cellInd, mixInd);
end

function highlightCell(cellInd, mixInd)
x = snrDb(cellInd,:);
y = r2(cellInd,:);
cols = interp1([1 0],[140 30 0;30 140 0]/140,linspace(0,1,11));
if (size(x,2)==16) % conc series included
  cols = [cols; [30 30 0; 60 60 0; 80 80 0; 100 100 0; 140 140 0]/140];
end
set(0,'CurrentFigure',clickMap);
if (isempty(hHighlightedCell) || ~ishandle(hHighlightedCell(1)))
  sfigure(clickMap);
  hold on;
  for i = 1:size(r2,2)
    hHighlightedCell(i) = plot3(x(i),y(i),10,'o','Color',cols(i,:),'MarkerSize',10,'MarkerFaceColor',cols(i,:));
    set(hHighlightedCell(i),'hittest','off');
  end
  hHighlightedCellMix = plot3(x(mixInd), y(mixInd),10,'wo','MarkerSize',20);
  set(hHighlightedCellMix,'hittest','off');
else
  for i = 1:size(r2,2)
    set(hHighlightedCell(i),'xdata',x(i),'ydata',y(i));
  end
  set(hHighlightedCellMix, 'xdata', x(mixInd), 'ydata', y(mixInd));
end
end
end