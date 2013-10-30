function MakeFiguresForPaper(whichPanels, varargin)
% function MakeFiguresForPaper(whichPanels, varargin)
%
% whichPanels = 1 => Reconstruction of some example KCs.
% whichPanels = 2 => Showing stats for KCs that could be reconstructed well.

p = inputParser;
p.addOptional('dataDir', 'originalData');
p.parse(varargin{:});

figDir  = GetDataDirForFigure(7);
dataDir = p.Results.dataDir;
currDir = GetCurrentDirFromPathString(fileparts(mfilename('fullpath')));

%% Panel 1: Showing some sample KC reconstructions
if (any(whichPanels==1))
  % Make the reconstructions figure
  data    = load(fullfile(figDir, currDir, dataDir, 'data100_t2_1_to_3_1.mat'));
  weights = load(fullfile(figDir, currDir, dataDir, 'prunedWeightsAdaptiveLasso.mat'));

  U = data.U;
  V = data.V;

  [numSamples, numPns] = size(U);
  [numSampels, numKcs] = size(V);

  B = weights.Results(:,1:174,1)';
  b0= weights.Results(:,175,1)';

  Vrec = bsxfun(@plus, U*B, b0);

  whichKcs = [35 84 21];
  cc  = KcReconstructionComparisonMetrics(V(:,whichKcs),Vrec(:,whichKcs),'corr',@Identity);
  ssen = KcReconstructionComparisonMetrics(V(:,whichKcs),Vrec(:,whichKcs),'ssen',@Identity);
  numConn= sum(B(:,whichKcs)~=0,1);

  sfigure(FindFigureCreate('Figure 7G: KC response reconstructions')); clf; 
  set(gcf,'Color',[1 1 1],'Resize','off','NumberTitle','off');
  ResizeFigure(gcf,10,8,'inches');
  Q = ComputeSubplotPositions(numel(whichKcs),1,[],0.05,0.02,0.00,0.025,0.075,0.01);
  set(gcf,'MenuBar','none','Toolbar','none');
  for i = 1:size(Q,1)
    subplotp(Q,i);
    plot(V(:,whichKcs(i)),'b','LineWidth',0.5);
    hold on;
    plot(Vrec(:,whichKcs(i)),'r','LineWidth',0.5);
    xlim([1 440]);
    ylim([-0.25 1.75]);
    set(gca,'xtick',([1 4 13 25 33 39 44]-1)*10+1,'ytick',-0.5:0.5:1.75);
    if (i~=size(Q,1))
      set(gca,'xticklabel',[]);
    else
      set(gca,'xticklabel',{'*high','1c','2c','3c','4c','5c','8c'});
    end
    set(gca,'ytick',0:0.5:1.75,'yticklabel',{'0','','1',''});
    grid off;
    set(gca,'FontSize',12);
    
    statsText = sprintf('#PNs = %d, SSE/SST = %1.2f, cc = %1.2f', numConn(i), ssen(i), cc(i));
    h = text(10,1.5,statsText,'FontSize',12);
  end
  legend('KC response', 'Rec. using PNs');
  xlabel('Odor',   'FontSize', 14);
  ylabel('Spikes', 'FontSize', 14);
end

%% Panel 2: Showing stats for the KCs that could be reconstructed well
if (any(whichPanels==2))
  % Now plot everything
  data = load(fullfile(figDir, currDir, dataDir, 'dataForKcRec.mat'));
  
  indSig      = data.indSig;
  ssenRec     = data.ssenRec;
  ssenRecPnSh = data.ssenRecPnSh;
  ssenRecKcSh = data.ssenRecKcSh;
  ccRec       = data.ccRec;
  ccRecPnSh   = data.ccRecPnSh;
  ccRecKcSh   = data.ccRecKcSh;

  sfigure(FindFigureCreate('Figure 7H: KC reconstructions statistics')); 
  clf; set(gcf,'Color',[1 1 1],'Resize','off','MenuBar','none','toolbar','none','NumberTitle','off');
  ResizeFigure(gcf,8,8,'inches');

  Q = ComputeSubplotPositions(1,2,[],0.075,0.075,0.02,0.075,0.02,0);
  subplotp(Q,1);
  % We'll plot the SSEn for 
  % (a) the well-reconstructed KCs,
  % (b) the well-reconstructed KCs reconstructed using shuffled PNs, pooled over shuffles
  % (c) all shuffled KCs reconstructed using un-shuffled PNs, pooled over KCs and shuffles.
  ssenData = {ssenRec(indSig), Columnize(ssenRecPnSh(:,indSig)), ssenRecKcSh(:)};
  muSsen = cellfun(@median, ssenData);
  lbSsen = cellfun(@(x) prctile(x,5),ssenData);
  ubSsen = cellfun(@(x) prctile(x,95),ssenData);
  colors = [name2rgb('red');name2rgb('darkolivegreen3');name2rgb('gray')];
  lgnd = [];
  for i = 1:3
    col = colors(i,:);
    lgnd(i) = line([i i],[lbSsen(i) ubSsen(i)],'LineWidth',1,'Color',col); hold on;
    plot(i,muSsen(i),'o','Color',col,'MarkerSize',5,'MarkerFaceColor',col,'MarkerEdgeColor',col); 
  end
  xlim([0 4]);
  ylim([0 1]);
  box on;
  set(gca,'xtick',[],'xticklabel',[],'ytick',[0:0.25:1],'yticklabel',{'0','','0.5','','1'},'FontSize',12);
  % Make a custom legend, because the default one is too big
  lgndY0 = 0.15;
  lgndDy = 0.05;
  lgndX0 = 0.2;
  lgndY = lgndY0;
  lgndText={'PNs, KCs','Shuff. PNs, KCs','PNs, Shuff. KCs'};
  for i = 1:3
    line([lgndX0 lgndX0+0.5],[lgndY lgndY],'Color',colors(i,:));
    text(lgndX0+0.6,lgndY-lgndDy*0.4,lgndText{i},'FontSize',6,'VerticalAlignment','Bottom');
    lgndY = lgndY - lgndDy;
  end
  title('SSE/SST','FontSize',16);
  
  % Plot the Correlation Coefficient Data
  subplotp(Q,2);
  ccData = {ccRec(indSig),  Columnize(ccRecPnSh(:,indSig)), ccRecKcSh(:)};
  muCc = cellfun(@nanmedian, ccData);
  lbCc = cellfun(@(x) prctile(x,5),ccData);
  ubCc = cellfun(@(x) prctile(x,95),ccData);

  for i = 1:3
    col = colors(i,:);
    line([i i],[lbCc(i) ubCc(i)],'LineWidth',1,'Color',col); hold on;
    plot(i,muCc(i),'o','Color',col,'MarkerSize',5,'MarkerFaceColor',col,'MarkerEdgeColor',col); 
  end
  set(gca,'yaxislocation','right');
  xlim([0 4]);
  ylim([0 1]);
  box on;
  set(gca,'xtick',[],'xticklabel',[],'ytick',[0:0.25:1],'yticklabel',{'0','','0.5','','1'},'FontSize',12);
  title('Corr. Coef.','FontSize',16);
  
end
