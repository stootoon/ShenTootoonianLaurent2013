function MakeFiguresForPaper(whichPanels, varargin)
% function MakeFiguresForPaper(whichPanels, ['dataDir' = 'originalData'])
%
% whichPanels=1 => Figure 5E: Transient proximity

p = inputParser;
p.addOptional('dataDir', 'originalData');
p.parse(varargin{:});

dataDir  = p.Results.dataDir;

figDir   = GetDataDirForFigure(5);
currDir  = GetCurrentDirFromPathString(fileparts(mfilename('fullpath')));
dataPath = fullfile(figDir, currDir, dataDir, 'Mfig_noduplicates');
data     = load(dataPath);
M        = data.M;

if (any(whichPanels==1))
  numBins = size(M.D,2);
  t = M.firstBinTime + (0:numBins-1)*M.binSize;
  whichCmp = 1;
  D = M.D(:,:,whichCmp);
  mu = mean(D,1);
  sd = std(D,[]);
  sem= sd/sqrt(size(D,1));
  
  sfigure(FindFigureCreate('Figure 5E: Transient proximity')); clf; 
  set(gcf,'Color',[1 1 1],'NumberTitle','off');
  ResizeFigure(gcf,8,5.33,'inches');
  Q = ComputeSubplotPosition(0.1,0.01,0.01,0.1);
  subplotp(Q,1);
  hp = patch([2.1 2.6 2.6 2.1],[0 0 1.2 1.2],'k');
  set(hp,'FaceAlpha',0.2,'EdgeColor','none');
  hl = PlotMeanWithFilledCiBand(t, mu, mu+sem, mu-sem, 'b', 1, 'b', 0.2);
  axis tight;
  xlim([1.5 5.0]);
  ylim([0 1.2]);
  set(gca,'xtick',1.5:0.5:5,'ytick',0:0.2:1.2,'FontSize',12);
  xlabel('Time(s)','FontSize',14);
  ylabel('Trajectory correlation distance','FontSize',14);
  box on;
end