% function MakeFiguresForPaper(whichFigures)
%
% Plots the Projection-Angle-Fraction related figures for the paper.
% WHICHFIGURES should be a vector of integers indicating the desired
% figures. The figures are:
% 
% 1: PAF wrt Citral (Fig. 3G)
% 2: PMF (Fig. S3A)
% 3: PAF wrt Octanol (Fig. S3B).

function MakeFiguresForPaper(whichFigures, varargin)
global Mfig numBs numJk numBins t plotMode

p = inputParser;
p.addOptional('dataDir', 'originalData');
p.parse(varargin{:});

dataDir = p.Results.dataDir;
thisDir = GetCurrentDirFromPathString(fileparts(mfilename('fullpath')));
figDir  = GetDataDirForFigure(3);

figPanels = {'3G', 'S3A', 'S3B'};

Mfig = LoadVarFromMatFileByName(fullfile(figDir, thisDir, dataDir, 'Mfig.mat'),'Mfig'); % Use the non-bootstrapped data.

if (isempty(Mfig.Ibs))
  plotMode  = 'mean_sem';
  plotDescr = 'Mean +/- SEM';
else
  plotMode  = 'mean_ci';
  plotDescr = 'Bootstrap bias-corrected and accelerated 5-95% CI';
end
  
numBs   = size(Mfig.Ibs,2)-1;
numBins = size(Mfig.PAFpbCit,2);
numJk   = size(Mfig.PAFpbCit,3)-numBs-1;
t       = Mfig.t0+(0:numBins-1)*Mfig.binSize;

if (any(whichFigures==1))
  % Plot the PAF timecourse figure
  figId = 1;
  figName = ['Figure ' figPanels{figId} ': PAF w.r.t. citral time course'];
  X = Mfig.PAFpbCit([1 6 11],:,:);
  doPlot(figId, figName, X, [name2rgb('ForestGreen'); name2rgb('Gold'); name2rgb('Red')]);
  disp(sprintf('%s for the PAF to citral for mostly citral ([30 140], green), 1:1 ([140 140], yellow), mostly octanol ([140 30], red). Xtick = 1.5:0.5:5, Ytick = [0:0.25:1].', plotDescr));  
end

if (any(whichFigures==2))
  % Plot the PMF timecourse SI figure
  figId = 2;
  figName = ['Figure ' figPanels{figId} ': PMF time course'];
  X = Mfig.PMFpb(6,:,:);
  lgnd = {'cit140:oct140'};
  doPlot(figId, figName, X, name2rgb('Gold'), lgnd);
  set(gca,'ytick',[0:0.25:1]);
  disp(sprintf('%s for the PMF for 1:1 mixture ([140 140], blue). Xtick = 1.5:0.5:5, Ytick = [0:0.25:1].', plotDescr));    
end

if (any(whichFigures==3))
  % Plot the PAF wrt Oct timecourse figure
  figId = 3;
  figName = ['Figure ' figPanels{figId} ': PAF w.r.t. octanol time course'];
  X = Mfig.PAFpbOct([1 6 11],:,:);
  doPlot(figId, figName, X, [name2rgb('Red'); name2rgb('Gold'); name2rgb('ForestGreen')]);
  disp(sprintf('%s for the PAF wrt octanol for mostly citral ([30 140], green), 1:1 ([140 140], yellow), mostly octanol ([140 30], red). Xtick = 1.5:0.5:5, Ytick = [0:0.25:1].', plotDescr));  
end

function doPlot(figId, figName, X, cols, lgnd)
global Mfig numBs numJk numBins t plotMode

if (nargin == 4)
  lgnd = {'cit140:oct30', 'cit140:oct140', 'cit30:oct140'};
end

sfigure(FindFigureCreate(figName)); clf; set(gcf,'Color',[1 1 1],'NumberTitle','off');
ResizeFigure(gcf,2.3*4,2*4,'inches');
Q = ComputeSubplotPosition(0.075,0.01,0.01,0.075);
subplotp(Q,1);

hp = patch([2.1 2.4 2.4 2.1],[0.005 0.005 0.995 0.995],name2rgb('lavender'),'EdgeColor','none');
Y = zeros(size(X,1), size(X,2), 3);
alpha = 0.05;
for i = 1:size(X,1) % numMixtureConds
  for j = 1:size(X,2) % numBins
    if (isequal(plotMode,'mean_sem'))
      mu = mean(X(i,j,:),3);
      sd = std(X(i,j,:));
      sem = sd/sqrt(size(X,3));
      Y(i,j,1) = mu - sem;
      Y(i,j,3) = mu + sem;
      Y(i,j,2) = mu;
    else
      uobs = squeeze(X(i,j,1));
      ubs  = squeeze(X(i,j,2:numBs+1));
      ujk  = squeeze(X(i,j,numBs+2:end));
      ci = ComputeBCaBootstrapConfidenceInterval(uobs, ubs, ujk, alpha);
      Y(i,j,1) = ci(1);
      Y(i,j,3) = ci(2);
      Y(i,j,2) = uobs;
    end
  end
end

hl = [];
for i = 1:size(Y,1)
  Z = squeeze(Y(i,:,:));
  hl(i) = PlotMeanWithFilledCiBand(t, Z(:,2), Z(:,1), Z(:,3), cols(i,:), 1, cols(i,:),0.2);    
end
set(legend(hl, lgnd{:}),'FontSize',14,'box','off');
set(gca,'xlim',[1.5 5], 'ylim',[0 1], 'ytick',[0:0.25:1],'xtick',0:0.5:5);
set(gca,'xticklabel', arrayfun(@num2str, get(gca,'xtick')-2, 'UniformOutput', false));
set(gca,'fontsize',12);
xlabel('Time (s)','FontSize', 14);
ylabel('PAF w.r.t. Citral','FontSize',14);

