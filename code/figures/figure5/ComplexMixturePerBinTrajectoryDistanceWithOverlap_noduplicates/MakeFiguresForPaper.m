function MakeFiguresForPaper(whichPanels, varargin)
% function MakeFiguresForPaper(whichPanels, ['dataDir' = 'originalData'])
%
% whichPanels=1 => Figure 5F:  Per-bin odor/trajectory relationship
% whichPanels=2 => Figure S5C: Per-bin odor/trajectory relationship
% whichPanels=3 => Figure S5D: Per-bin odor/trajectory data

p = inputParser;
p.addOptional('dataDir', 'originalData');
p.parse(varargin{:});

dataDir = p.Results.dataDir;

figDir   = GetDataDirForFigure(5);
currDir  = GetCurrentDirFromPathString(fileparts(mfilename('fullpath')));
dataPath = fullfile(figDir, currDir, dataDir, 'Mfig_noduplicates.mat');
data     = load(dataPath);
M        = data.M;

if (any(whichPanels==1))  
  % Plot the time course of the spearman correlation in blue
  % and the shuffled time course in red.

  numBins    = size(M.sp,4);
  numSamples = size(M.sp,3);
  t          = M.firstBinTime + (0:numBins-1)*M.binSize;

  sfigure(FindFigureCreate('Figure 5F: Per-bin odor/trajectory relationship')); clf; 
  set(gcf,'Color',[1 1 1],'NumberTitle','off');
  ResizeFigure(gcf,8,5.33,'inches');
  
  Q = ComputeSubplotPosition(0.1,0.01,0.01,0.1);
  subplotp(Q,1);
  
  whichTrajDist = 1; % correlation dist
  whichOdorDist = 2; % intersection dist

  mu  = squeeze(mean(M.sp(whichTrajDist, whichOdorDist, :, :),3));
  sd  = squeeze( std(M.sp(whichTrajDist, whichOdorDist, :, :),[], 3));
  sem = sd/sqrt(numSamples);

  mush  = squeeze(mean(M.spsh(whichTrajDist, whichOdorDist, :, :),3));
  sdsh  = squeeze( std(M.spsh(whichTrajDist, whichOdorDist, :, :),[], 3));
  semsh = sd/sqrt(numSamples);

  muPnSh  = squeeze(mean(M.spPnSh(whichTrajDist, whichOdorDist, :, :, 1),3));
  sdPnSh  = squeeze( std(M.spPnSh(whichTrajDist, whichOdorDist, :, :, 1),[], 3));
  semPnSh = sd/sqrt(numSamples);
  
  odorOnsetTime = 2;
  
  hp = patch([2.1 2.6 2.6 2.1] - odorOnsetTime,[-0.2 -0.2 1 1],name2rgb('black'));
  set(hp,'FaceAlpha',0.2,'EdgeColor','none'); hold on;
  hl     = PlotMeanWithFilledCiBand(t-odorOnsetTime, mu, mu+sem, mu-sem, 'b', 1, 'b', 0.2);
  hlsh   = PlotMeanWithFilledCiBand(t-odorOnsetTime, mush, mush+semsh, mush-semsh, 'r', 1, 'r', 0.2);  
  hlpnsh = PlotMeanWithFilledCiBand(t-odorOnsetTime, muPnSh, muPnSh+semPnSh, muPnSh-semPnSh, 'k', 1, 'k', 0.2);
  box on;  
  axis tight;
  ylim([-0.2 1]);
  set(gca,'xtick',[-0.5:0.5:3.0],'ytick',[-0.2:0.2:1],'FontSize',12);
  xlabel('Time (s)','FontSize', 14);
  ylabel('Spearman rank correlation', 'FontSize', 14);
  legend([hl hlsh hlpnsh], 'unperturbed','odor labels shuffled','PN identities shuffled');
end

if (any(whichPanels==2))  
  % Plot the time course of the spearman correlation in blue
  % and the shuffled time course in red, for each of the odor distance functions
  
  numBins = size(M.sp,4);
  numSamples = size(M.sp,3);
  t = M.firstBinTime + (0:numBins-1)*M.binSize;

  sfigure(FindFigureCreate('Figure S5C: Per-bin odor/trajectory relationship')); clf; 
  set(gcf, 'Color', [1 1 1], 'NumberTitle', 'off');
  ResizeFigure(gcf, 12, 5, 'inches');
  
  Q = ComputeSubplotPositions(1,2,[],0.075,0.01,0.02,0.1,0.12,0.0);
  
  whichTrajDist  = 1; % correlation dist
  whichOdorDists = [1 3];
  odorOnsetTime  = 2;
  ttls = {'Odor Braun-Blanquet distance', 'Odor cosine distance'};
  
  for i = 1:numel(whichOdorDists)
    subplotp(Q,i);
    whichOdorDist = whichOdorDists(i);
  
    mu  = squeeze(mean(M.sp(whichTrajDist, whichOdorDist, :, :),3));
    sd  = squeeze(std(M.sp(whichTrajDist, whichOdorDist, :, :),[], 3));
    sem = sd/sqrt(numSamples);

    mush  = squeeze(mean(M.spsh(whichTrajDist, whichOdorDist, :, :),3));
    sdsh  = squeeze(std(M.spsh(whichTrajDist, whichOdorDist, :, :),[], 3));
    semsh = sd/sqrt(numSamples);

    muPnSh  = squeeze(mean(M.spPnSh(whichTrajDist, whichOdorDist, :, :,1),3));
    sdPnSh  = squeeze(std(M.spPnSh(whichTrajDist, whichOdorDist, :, :,1),[], 3));
    semPnSh = sd/sqrt(numSamples);

    hp = patch([2.1 2.6 2.6 2.1] - odorOnsetTime, [-0.2 -0.2 1 1], name2rgb('black'));
    set(hp,'FaceAlpha',0.2,'EdgeColor','none'); hold on;
    hl     = PlotMeanWithFilledCiBand(t - odorOnsetTime, mu, mu+sem, mu-sem, 'b', 1, 'b', 0.2);
    hlsh   = PlotMeanWithFilledCiBand(t - odorOnsetTime, mush, mush+semsh, mush-semsh, 'r', 1, 'r', 0.2);
    hlpnsh = PlotMeanWithFilledCiBand(t - odorOnsetTime, muPnSh, muPnSh+semPnSh, muPnSh-semPnSh, 'k', 1, 'k', 0.2);
    box on;  
    axis tight;
    ylim([-0.2 1]);
    set(gca,'xtick',[1:0.5:10]-odorOnsetTime,'ytick',[-0.2:0.2:1],'FontSize',12);
    if (i==1)
      xlabel('Time(s)', 'FontSize', 14);
      ylabel('Spearman rank correlation','FontSize', 14);
      legend([hl hlsh hlpnsh], 'unperturbed','odor labels shuffled','PN identities shuffled');
    else
      set(gca,'xticklabel',[],'yticklabel',[]);
    end
    title(ttls{i},'FontSize',16);
  end
end

if (any(whichPanels==3))  
  % Plot the data used in the spearman correlation computations at various time points.

  numBins    = size(M.sp,4);
  t          = M.firstBinTime + (0:numBins-1)*M.binSize;

  whichTimes = [1.9 2.2 2.9 4.0];
  whichBins  = arrayfun(@(wt) find(abs(t-wt)<0.01), whichTimes);
  
  sfigure(FindFigureCreate('Figure S5D: Per-bin odor/trajectory data')); clf; 
  set(gcf,'Color',[1 1 1],'NumberTitle','off');
  ResizeFigure(gcf,12,4,'inches');
  
  Q = ComputeSubplotPositions(1,4,[],0.07,0.02,0.01,0.025,0.15,0.0);
  
  whichTrajDist = 1; % Correlation distance
  whichOdorDist = 2; % Intersection distance
  
  dodors = M.dodors(whichOdorDist, :);
  dx     = min(diff(unique(dodors)));

  dr = reshape(M.dr(whichTrajDist, whichOdorDist, :, whichBins,:), size(M.dr,3), numel(whichBins), size(M.dr,5));

  x = dodors;
  odorOnsetTime = 2;
  for i = 1:numel(whichBins)
    subplotp(Q,i);
    Y = squeeze(dr(:,i,:));
    cols = jet(7);
    for j = 1:size(Y,1)
      plot(x+randn(size(x))*dx/60+j*dx/10-dx/5, Y(j,:),'.','Color',cols(j,:));
      hold on;
    end
    xlim([0 1.2]);
    ylim([0 1.2]);
    set(gca,'xtick',0:0.2:1.2,'ytick',0:0.2:1.2,'FontSize',12);
    if (i>1)
      set(gca,'xticklabel',[],'yticklabel',[]);
    else
      xlabel('Odor Jaccard distance','FontSize',14);
      ylabel('Trajectory correlation distance','FontSize',14);
    end
    text(0.6,0.05,sprintf('T = %1.1f', t(whichBins(i))-odorOnsetTime),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize', 16);
  end
end
    
