function MakeFiguresForPaper(whichPanels, varargin)
% function MakeFiguresForPaper(whichPanels, ['dataDir'='originalData'])
%
% whichPanels=1 => Figure 5D:  Global odor/trajectory relationship
% whichPanels=2 => Figure S5A: Global odor/trajectory data
% whichPanels=3 => Figure S5B: Global odor/trajectory relationships (different odor metrics)

p = inputParser;
p.addOptional('dataDir', 'originalData');
p.parse(varargin{:});

dataDir = p.Results.dataDir;

figDir   = GetDataDirForFigure(5);
thisDir  = GetCurrentDirFromPathString(fileparts(mfilename('fullpath')));
dataPath = fullfile(figDir, thisDir, dataDir, 'Mfig_noduplicates.mat');
data     = load(dataPath);
M        = data.M;

if (any(whichPanels==1))  
  % Plot a bar graph showing the mean and SEM of the spearman rank correlations.
  sfigure(FindFigureCreate('Figure 5D: Global odor/trajectory relationship')); clf;
  set(gcf,'Color',[1 1 1],'NumberTitle','off');
  ResizeFigure(gcf,8,8,'inches');
  Q = ComputeSubplotPosition(0.1,0.01,0.01,0.01);
  subplotp(Q,1);

  whichTrajDist = 1;
  whichOdorDist = 2;

  mu = [mean(M.spsh(whichTrajDist, whichOdorDist,:),3);
        mean(M.spbl(whichTrajDist, whichOdorDist,:),3);
        mean(M.sp(whichTrajDist, whichOdorDist,:),3)]';
  
  sd = [std(M.spsh(whichTrajDist,   whichOdorDist,:), [], 3);
        std(M.spbl(whichTrajDist, whichOdorDist,:), [], 3);
        std(M.sp(whichTrajDist, whichOdorDist,:), [], 3)]';
  
  sem = sd/sqrt(size(M.sp,3));
  
  hb = BarPlot(mu, mu-sem, mu+sem, 0.4, [0.5 0.5 0.5; 1 0 0; 0 0 1;], 1);
  ylim([-0.2 1]);
  set(gca,'xtick',[],'FontSize',12);
  xlim([0 4]);
  h = line([0 4],[0 0]);
  set(h,'LineStyle','--','Color','k');
  ylabel('Spearman rank correlation','FontSize',14);
  legend(hb(end:-1:1), 'response window','baseline','response window (odor labels shuffled','Location','NorthWest');
  box on;
end

  
if (any(whichPanels==2))
  % Plot the data points for each trial
  sfigure(FindFigureCreate('Figure S5A: Global odor/trajectory data')); clf; 
  set(gcf,'Color',[1 1 1],'NumberTitle','off');
  ResizeFigure(gcf,8,8,'inches');
  Q = ComputeSubplotPosition(0.1,0.01,0.01,0.1);
  subplotp(Q,1);
  whichTrajDist = 1; % Correlation distance
  whichOdorDist = 2; % Intersection distance
  
  dodors = Columnize(M.dodors(whichTrajDist, whichOdorDist, 1, :));
  Dresp  = squeeze(M.dresp(whichTrajDist, whichOdorDist, :,:));
  
  dx = min(diff(unique(dodors)));
  
  cols = jet(7);
  for i = 1:size(Dresp,1)
    plot(dodors + randn(size(dodors))*dx/30+i*(dx/10)-dx/5, Dresp(i,:),'.','Color',cols(i,:),'MarkerSize',5); hold on;
  end

  xlim([0 1.2]);
  ylim([0 1.2]);
  set(gca,'xtick',0:0.2:1.2,'ytick',0:0.2:1.2,'FontSize', 12);
  xlabel('Odor Jaccard distance','FontSize', 14);
  ylabel('Trajectory correlation distance', 'FontSize', 14);
  TightenAxesToFigure;
end


if (any(whichPanels == 3))
  % Plot a bar graph for each of the distance function pairs.
  numTrajDists = 1; % Only using correlation distance
  whichOdorDists = [1 3]; % Use only the overlap and cosine distances, since the intersection distance is already in the main figure.
  numOdorDists = numel(whichOdorDists);

  sfigure(FindFigureCreate('Figure S5B: Global odor/trajectory relationships (different odor metrics)')); clf; 
  set(gcf,'Color',[1 1 1],'NumberTitle','off');
  ResizeFigure(gcf,8,4,'inches');
  Q = ComputeSubplotPositions(numTrajDists, numOdorDists,[],0.1,0.01,0.01,0.1,0.05,0.0);

  ttls = {'Odor Braun-Blanquet distance','Odor cosine distance'};
  for i = 1:1 % Only use correlation distance for trajectories. 
    for j = 1:numOdorDists
      subplotp(Q,j+(i-1)*numOdorDists);
      
      mu = [mean(M.spsh(i, whichOdorDists(j),:), 3);
            mean(M.spbl(i, whichOdorDists(j),:), 3);
            mean(M.sp(i, whichOdorDists(j),:),   3)]';
      
      sd = [std(M.spsh(i, whichOdorDists(j),:), [], 3);
            std(M.spbl(i, whichOdorDists(j),:), [], 3);
            std(M.sp(i, whichOdorDists(j),:),   [], 3)]';
  
      sem = sd/sqrt(size(M.sp,3));
  
      BarPlot(mu, mu-sem, mu+sem, 0.4, [0.5 0.5 0.5; 1 0 0; 0 0 1;], 1);
      if (j>1)
        set(gca,'yticklabel',[]);
      else
        ylabel('Spearman rank correlation', 'FontSize', 14);
      end
      ylim([-0.2 1]);
      set(gca,'xtick',[],'FontSize',12);
      xlim([0 4]);
      h = line([0 4],[0 0]);
      set(h,'LineStyle','--','Color','k');
      box on;
      title(ttls{j},'FontSize',15);
    end
  end
end
