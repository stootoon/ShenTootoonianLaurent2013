function MakeFiguresForPaper(varargin)
% function MakeFiguresForPaper(varargin)

p = inputParser;
p.addOptional('dataDir','originalData');
p.parse(varargin{:});
dataDir = p.Results.dataDir;
currDir = GetCurrentDirFromPathString(fileparts(mfilename('fullpath')));
figDir  = GetDataDirForFigure(7);

sfigure(FindFigureCreate('Figure S7A: Time course of PN and KC responses at each mixuture level.')); clf;
set(gcf,'Color',[1 1 1],'Resize','off','NumberTitle','off');
ResizeFigure(gcf,12,7.5,'inches');

Q = ComputeSubplotPositions(6,7,[],0.1,0.01,0.01,0.05,0.075,0.05);

mlPlotOrder = [1 0 2:5 8];
Data = LoadVarFromMatFileByName(fullfile(figDir, currDir, dataDir, 'dataPnKcResponseTimeCourse.mat'),'Data');

%% First plot the Mfrvsmls
disp('Plotting Mfrvsml...');
pn1CmpVal = 0;
kc1CmpVal = 0;
ttls = {'1-','4x 1-', '2-', '3-', '4-', '5-', '8-'};
for i = 1:numel(mlPlotOrder)
  indMl = find(Data.mixtureLevels == mlPlotOrder(i));
  pnData = Data.pnMfrvsml{indMl};
  kcData = Data.kcMfrvsml{indMl};

  if (mlPlotOrder(i)==1)
    pn1CmpVal = max(mean(pnData,2));
    kc1CmpVal = max(mean(kcData,2));    
  end
  
  numBins = Data.pnSz50(2);
  binSize = 0.05;
  t = (0:numBins-1)*binSize;
  xl = [0 3];
  xtick = 0:3;
  
  subplotp(Q,i);
  yl = [1.5 4.5];
  ytick = yl;
  doSinglePlot(pnData/binSize, t, xl, yl, xtick, ytick, [0, i==1])
  hold on;
  HorizontalLine(gca, pn1CmpVal/binSize, 'Color', name2rgb('gray70'));
  if (i==1)
    ylabel('PNs','FontSize',14);
    set(gca,'FontSize',12);
  end
  title(ttls{i},'FontSize',16);
  
  subplotp(Q,i+7);
  yl = [0 1];
  ytick = yl;
  doSinglePlot(kcData/binSize, t, xl, yl, xtick, ytick, [0, i==1])  
  hold on;
  HorizontalLine(gca, kc1CmpVal/binSize, 'Color', name2rgb('gray70'));
  if (i==1)
    h = ylabel('KCs','FontSize',14);
    c = copyobj(h,get(h,'Parent'));
    pos = get(c,'Position'); pos(1:2) = [-1.5 0];
    set(c,'Position',pos,'HorizontalAlignment','left','string', 'Mean Firing Rate (Hz)','FontSize',16);
    set(gca,'FontSize',12);
  end

end

%% Now plot the silent cells
disp('Plotting Silence Percentage...');
for i = 1:numel(mlPlotOrder)
  indMl = find(Data.mixtureLevels == mlPlotOrder(i));
  pnData = Data.pnSilent{indMl};
  kcData = Data.kcSilent{indMl};

  if (mlPlotOrder(i)==1)
    pn1CmpVal = max(mean(pnData,2));
    kc1CmpVal = min(mean(kcData,2));
  end

  numBins = Data.pnSz100(2);
  binSize = 0.100;
  t = (0:numBins-1)*binSize;

  xl = [0 3];
  xtick = 0:3;
  
  subplotp(Q,i+2*7);
  yl = [30 90];
  ytick = [30 60 90];
  doSinglePlot(pnData*100, t, xl, yl, xtick, ytick, [0, i==1]);
  hold on;
  HorizontalLine(gca, pn1CmpVal*100, 'Color', name2rgb('gray70'));
  if (i==1)
    ylabel('PNs','FontSize',14);
    set(gca,'FontSize',12);
  end
  
  subplotp(Q,i+3*7);
  yl = [60 100];
  ytick = yl;
  doSinglePlot(kcData*100, t, xl, yl, xtick, ytick, [0, i==1]);
  hold on;
  HorizontalLine(gca, kc1CmpVal*100, 'Color', name2rgb('gray70'));
  if (i==1)
    h = ylabel('KCs','FontSize',14);
    c = copyobj(h,get(h,'Parent'));
    pos = get(c,'Position'); pos(1:2) = [-1.5 80];
    set(c,'Position',pos,'HorizontalAlignment','left','string', 'Silent Cells (%)','FontSize',16);
    set(gca,'FontSize',12);
  end
end

%% Now plot the responsive cells
disp('Plotting Responsive Percentage...');
for i = 1:numel(mlPlotOrder)
  indMl = find(Data.mixtureLevels == mlPlotOrder(i));
  pnData = Data.pnResponsive{indMl};
  kcData = Data.kcResponsive{indMl};
 
  numBins = Data.pnSz50(2);
  binSize = 0.050;
  t = (0:numBins-1)*binSize;

  xl = [0 3];
  xtick = 0:3;
  
  subplotp(Q,i+4*7);
  yl = [0 14];
  ytick = [0 7 14];
  doSinglePlot(pnData*100, t, xl, yl, xtick, ytick, [0, i==1]);
  text(numBins*2/3*binSize, 2/3*yl(2), sprintf('%d %%', round(Data.pnCumResp(indMl)*100)),'FontSize', 12);
  if (i==1)
    ylabel('PNs','FontSize',14);
    set(gca,'FontSize',12);
  end
  
  subplotp(Q,i+5*7);
  yl = [0 1.2];
  ytick = [0 0.6 1.2];
  doSinglePlot(kcData*100, t, xl, yl, xtick, ytick, [i==1, i==1]);
  text(numBins*2/3*binSize, 2/3*yl(2), sprintf('%1.1f %%', Data.kcCumResp(indMl)*100), 'FontSize', 12);
  if (i==1)
    h = ylabel('KCs','FontSize',14);
    c = copyobj(h,get(h,'Parent'));
    pos = get(c,'Position'); pos(1:2) = [-1.5 0];
    set(c,'Position',pos,'HorizontalAlignment','left','string', 'Responsive Cells (%)','FontSize',16);
    set(gca,'FontSize',12);
    xlabel('Time (s)', 'FontSize', 14);
  end
end

function doSinglePlot(data, t, xl, yl, xtick, ytick, makeLabels)
ghostColor = name2rgb('gray90');
traceColor = name2rgb('black');

hp = patch([0 0.5 0.5 0]+0.1,[yl(1) yl(1) yl(2) yl(2)],name2rgb('magenta'));
set(hp,'EdgeColor','none','FaceAlpha',0.1);
hold on;
plot(t, data,'Color',ghostColor);
if (~IsVector(data))
  mu = mean(data,2);
else
  mu = data; 
end
hold on;
plot(t,mu,'Color',traceColor);
ylim(yl);
xlim(xl);
set(gca,'xtick',xtick,'xticklabel',[],'ytick',ytick,'yticklabel',[]);

if (makeLabels(1))
  set(gca,'xticklabel',arrayfun(@num2str, xtick, 'UniformOutput', false));
end

if (makeLabels(2))
  set(gca,'yticklabel',arrayfun(@num2str, ytick, 'UniformOutput', false));
end
