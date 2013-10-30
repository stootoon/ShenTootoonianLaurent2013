function MakeFiguresForPaper(whichFigures, varargin)
% MakeFiguresForPaper(whichFigures)
%
% Makes the figures for the trajectory morphs. 
% 
% Figure 1: Log probabilities of the models relative to the linear
% model for the fit to the log of the concentration ratio vs
% normalized correlation distance to citral.
%
% Figure 2: The SI for figure 1, but now using percentage of the
% mixture rather log of the concentration ratio.
%
% Figure 3: Another SI, this time showing the data and the fits for
% one trial of figure 1,

% Plot the mean and sem of the logps for each of the models of distance relative to citral
p = inputParser;
p.addOptional('dataDir','originalData');
p.parse(varargin{:});

dataDir  = p.Results.dataDir;
currDir  = GetCurrentDirFromPathString(fileparts(mfilename('fullpath')));
figDir   = GetDataDirForFigure(3);

addpath(fileparts(mfilename('fullpath')));

if (any(whichFigures==1))
  figId = 1;
  figName = 'LogProbsRelLinearForLogConcRatioVsNormalizedCorrDistFromCitral';
  figTitle= 'Figure 3I: Trajectory evolution models';
  load(fullfile(figDir, currDir, dataDir, 'Mfiglog'));
  doPlotTimecourse(figId, figTitle, Mfig, 'correlation');
end

if (any(whichFigures==2))
  figId = 101;
  figName = 'LogProbsRelLinearForPcConcVsCorrDistFromCitral';
  figTitle= 'Figure S3C: Trajectory evolution models (i.v. = frac. oct.)';
  load(fullfile(figDir, currDir, dataDir, 'Mfigpc'));
  doPlotTimecourse(figId, figTitle, Mfig, 'correlation');
end

if (any(whichFigures == 3))
  figId = 102;
  figName = 'FitsAtDifferentTimePointsForLogConcRatioVsNormalizedCorrDistFromCitral';
  figTitle = 'Figure S3D: Model fits at different time points (i.v. Log10([Oct]/[Cit]))'; 
  load(fullfile(figDir, currDir, dataDir, 'Mfiglog'));
  Results = Mfig.Results(5);
  x  = Results.x;
  Y  = Results.Y;
  numBins = size(Y,2);
  yr = Results.yr;
  t  = Mfig.t0 + (0:numBins-1)*Mfig.binSize;
  whichTimes = [1.9 2.2 2.3 2.5];
  whichBins = arrayfun(@(tt) find(t==tt,1), whichTimes);
  whichTrial = 4;
  Y = squeeze(Y(whichTrial, whichBins, :));
  xl = [-1:0.25:1];
  hp = doPlotFitAtTimePoint(figId, figTitle, x, Y, yr, xl); 
  xlabel(hp(1),'Log_1_0([Oct]/[Cit])','FontSize',14);
  ylabel(hp(1),'Correlation Distance','FontSize',14);
  set(hp(1), 'xticklabel', arrayfun(@(x) when(abs(mod(x,0.5))<0.1, sprintf('%1.1f', x),'') , get(gca,'xtick'),'UniformOutput', false));
  set(hp(1), 'yticklabel', arrayfun(@(x) when(abs(mod(x,0.5))<0.1, sprintf('%1.1f', x),'') , get(gca,'ytick'),'UniformOutput', false));
  arrayfun(@(i) title(hp(i), sprintf('T = %1.1f s', whichTimes(i) - 2),'FontSize',16), 1:numel(hp));
end

if (any(whichFigures == 4))
  figId = 103;
  figName  = 'FitsAtDifferentTimePointsForPcConcVsNormalizedCorrDistFromCitral';
  figTitle = 'Figure S3D: Model fits at different time points (i.v. Fraction Octanol)'; 
  load(fullfile(figDir, currDir, dataDir, 'Mfigpc'));
  Results = Mfig.Results(5);
  x  = Results.x;
  Y  = Results.Y;
  numBins = size(Y,2);
  yr = Results.yr;
  t  = Mfig.t0 + (0:numBins-1)*Mfig.binSize;
  whichTimes = [1.9 2.2 2.3 2.5];
  whichBins = arrayfun(@(tt) find(t==tt,1), whichTimes);
  whichTrial = 4;
  Y = squeeze(Y(whichTrial, whichBins, :));
  xl = [0:0.25:1];
  hp = doPlotFitAtTimePoint(figId, figTitle, x, Y, yr, xl);  
  xlabel(hp(1),'Fraction Octanol','FontSize',14);
  ylabel(hp(1),'Correlation Distance','FontSize',14);
  set(hp(1), 'xticklabel', arrayfun(@(x) when(abs(mod(x,0.25))<0.05, sprintf('%1.2f', x),'') , get(gca,'xtick'),'UniformOutput', false));
  set(hp(1), 'yticklabel', arrayfun(@(x) when(abs(mod(x,0.5))<0.1, sprintf('%1.1f', x),'') , get(gca,'ytick'),'UniformOutput', false));
  arrayfun(@(i) title(hp(i), sprintf('T = %1.1f s', whichTimes(i) - 2),'FontSize',16), 1:numel(hp));
end

function hp = doPlotFitAtTimePoint(figId, figName, x, Y, yr, xl)
xr = [min(x) max(x)];
nmc = 10000;
posteriorFuncs = {@(x,y,xr,yr) ComputePosteriorOfConstantModelDirectly(x,y,xr,yr),...
                  @(x,y,xr,yr) ComputePosteriorOfLinearModel(x,y,xr,yr),...
                  @(x,y,xr,yr) ComputePosteriorOfSingleStepModelDirectly(x,y,xr,yr),...
                  @(x,y,xr,yr) ComputePosteriorOfTwoStepModel(x,y,xr,yr,nmc)};
numPoints = size(Y,1);
logp = zeros(1,4);
D = {};
sfigure(FindFigureCreate(figName)); clf; set(gcf,'Color',[1 1 1]); 
ResizeFigure(gcf,12,3.5,'inches'); set(gcf,'NumberTitle','off');
Q = ComputeSubplotPositions(1,numPoints,[],0.05,0.01,0.01,0.1,0.15,0.0);
yl = [min(Y(:)) max(Y(:))];
xi = linspace(xl(1),xl(end),100);

for i = 1:numPoints
  hp(i) = subplotp(Q,i);
  plot(x,Y(i,:),'ko','MarkerFaceColor','k'); hold on;
  set(gca,'yticklabel', [],'xticklabel',[]);
  ylim([0 2]);
  xlim([xl(1) xl(end)]);
  set(gca,'ytick',[0:0.25:2]);
  set(gca,'xtick',xl);
  for j = 1:numel(posteriorFuncs)
    [logp(j), D{j}] = posteriorFuncs{j}(x,Y(i,:)',xr,yr);
  end
  logpText = sprintf('%1.2f  ', logp*log10(exp(1)));
  [foo, im] = sort(logp,'descend');
  bestModel = im(1);
  secondBestModel = im(2);
  plot(xi, D{secondBestModel}.fval(xi), 'Color',[0.75 0.75 0.75]);
  plot(xi, D{bestModel}.fval(xi), 'Color', [1 0 0]);
  text(xl(1),1.85,logpText,'FontSize',12);
end
disp(sprintf('Data at various time points and the best (red) and secondbest (gray) fits. Log10 probs of all models are indicated. Xticks are at [%1.1f:%1.1f:%1.1f] and yticks are at [0:0.25:2]. Log10 probs of all models are indicated at the top left.',xl(1), xl(2)-xl(1), xl(end)));

function doPlotTimecourse(figId, figName, Mfig, whichDistance)
numBins = size(Mfig.Results(end).logp,2);
t = Mfig.t0+(0:numBins-1)*Mfig.binSize;
switch(whichDistance)
 case 'correlation'
  ind = 5;
  ydataName = 'normalized correlation distance to citral';
 case 'euclidean'
  ind = 4;
  ydataName = 'normalized euclidean distance to citral.';
 case 'PAF'
  ind = 6;
  ydataName = 'PAF to citral';    
end

switch(Mfig.coords)
 case 'log'
  xdataName = 'log concentration ratio';
 case 'pc'
  xdataName = 'percent concentration of octanol';
end
  
X = (Mfig.Results(ind).logp)*log10(exp(1)); % X is in samples x bins x models format
X = bsxfun(@minus, X, X(:,:,2)); % Relative to linear model
if (isempty(Mfig.Ibs))
  Y = zeros(3, size(X,2), size(X,3));
  for i = 1:size(X,2)
    for j = 1:size(X,3)
      mu = mean(X(:,i,j));
      sd = std(X(:,i,j));
      sem = sd/sqrt(size(X,1));
      Y(2,i,j) = mu;
      Y(1,i,j) = mu - sem;
      Y(3,i,j) = mu + sem;
      end
  end
  plotMode = 'Mean +/- SEM';
else
  pc = [5 50 95];
  Y = prctile(X, pc);
  plotMode = 'Boostrap [5 - 95] Percentile Intervals';
end

sfigure(FindFigureCreate(figName)); clf;
ResizeFigure(gcf,8,8,'inches');
set(gcf,'Color',[1 1 1],'numberTitle','off');
Q = ComputeSubplotPosition(0.11,0.02,0.01,0.075);
subplotp(Q,1);

hp = patch([2.1 2.4 2.4 2.1],[-4.97 -4.97 2 2],name2rgb('lavender'), 'EdgeColor','none','FaceAlpha',1); hold on;
models = {'const', 'lin','1step','2step'};
% We're plotting the logs relative to the linear model. 
% We'll color the traces gray (constant model), blue (1step), red (2step)
cols = [name2rgb('gray75'); name2rgb('blue'); name2rgb('red'); name2rgb('orange');]; 
hl = [];
for i = 1:4
  Z = squeeze(Y(:,:,i));
  hl(i) = PlotMeanWithFilledCiBand(t, Z(2,:),Z(3,:),Z(1,:), cols(i,:),1,cols(i,:),0.2); 
  set(hl(i),'LineWidth',2);
end
xlim([1.5 3.5]);
set(gca,'xtick',1.5:0.5:3.5);
set(gca,'xticklabel',arrayfun(@(x) num2str(x-2), get(gca,'xtick'), 'UniformOutput', false));
ylim([-5 2]);
set(gca,'ytick',[-5:2]);
xlabel('Time (s)','FontSize',14);
ylabel({'Log_1_0 Prob(Model)','w.r.t. linear model'},'FontSize', 14);
set(gca,'FontSize',12);
set(legend(hl([2 1 3 4]),'linear model','constant','one-step','two-step','location','SouthEast'),'box','off','FontSize',14);
disp([plotMode ' of Log10 prob of const (gray), 1step (red), 2step (orange) relative to linear (blue), for fitting ' xdataName ' to ' ydataName '. Xticks are at [1.5:0.5:3.5] and yticks are at [-5:2]']);
