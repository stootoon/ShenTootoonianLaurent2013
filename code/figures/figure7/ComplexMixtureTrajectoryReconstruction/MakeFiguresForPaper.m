function MakeFiguresForPaper(whichPanels, varargin)
% MakeFiguresForPaper(whichPanels, ['dataDir' = 'originalData'])
%
% whichPanels = 1 => Figure 7J: PN trajectory reconstructions
% whichPanels = 2 => Figure 7K: PN trajectory reconstruction statistics
p = inputParser;
p.addOptional('dataDir', 'originalData');
p.parse(varargin{:});

figDir  = GetDataDirForFigure(7);
dataDir = p.Results.dataDir;
currDir = GetCurrentDirFromPathString(fileparts(mfilename('fullpath')));

numBinsPerOdor = 10;
numOdors = 44;
convert2cbo = @(X) permute(reshape(X,numBinsPerOdor,numOdors,[]),[3 1 2]);

NEWFIG = @(i) eval(strrep(sprintf('sfigure(%d); clf; set(gcf,"Color",[1 1 1],"Resize","off"); ResizeFigure(gcf,8,8,"inches");', i),'"',''''));

if (any(whichPanels == 1))  % Projection onto the PCA axes
  % Showing the PCA reconstructions
  data = load(fullfile(figDir, currDir, dataDir, 'data100_t2_1_to_3_1.mat'));
  conn = load(fullfile(figDir, currDir, dataDir, 'prunedWeightsAdaptiveLasso.mat'));

  B = conn.Results(:,1:206,1)';
  B0= conn.Results(:,207,1)';
  U = data.U;
  V = data.V;
  [numSamples,numPns] = size(U);
  [numSamples,numKcs] = size(V);

  [Urec,Upn] = ReconstructPnTrajectories(V,U,B,B0);
  
  [pc,coefPn,latent] = princomp(U);
  whichOdor = 25;
  Zpn = Upn(:,:,whichOdor);
  Zrec= Urec(:,:,whichOdor);

  coefPn = pc'*Zpn;
  coefRec = pc'*Zrec;

  sfigure(FindFigureCreate('Figure 7J: PN trajectory reconstructions')); clf;
  set(gcf,'Color',[1 1 1],'Resize','off','NumberTitle','off'); ResizeFigure(gcf,9,7.5,'inches');

  lineWidth = 1;
  
  lgnd = [];
  lgnd(1) = plot3(coefPn(1,:),coefPn(2,:),coefPn(3,:),'bo-','LineWidth',lineWidth,'MarkerSize',3); hold on;
  plot3(coefPn(1,1),coefPn(2,1),coefPn(3,1),'o','Color','b','MarkerFaceColor','b','MarkerSize',5);
  text(coefPn(1,1), coefPn(2,1), coefPn(3,1), 'PN trajectory', 'Color','b','FontSize',12);
  
  lgnd(2) = plot3(coefRec(1,:),coefRec(2,:),coefRec(3,:),'ro-','LineWidth',lineWidth,'MarkerSize',3); 
  plot3(coefRec(1,1),coefRec(2,1),coefRec(3,1),'o','Color','r','MarkerFaceColor','r','MarkerSize',5);
  text(coefRec(1,1), coefRec(2,1), coefRec(3,1), 'Rec. using KCs', 'Color','r','FontSize',12);
  
  dataKcSh = load(fullfile(figDir, currDir, dataDir, 'lassoWeightsKcShuffle.mat'));
  [Vsh,VshCbo,Ush,UshCbo,Bsh,b0sh] = CollectDataForShuffledKcs(dataKcSh.Results);  
  numShuffles = size(Vsh,3);

  numShufflesToPlot = 3;
  % Find the shuffles which have the highest correlation for this odor.
  dc = zeros(1,numShuffles);
  for i = 1:numShuffles
    dc(i) = min(corrdist(Urec(:,:,whichOdor),UshCbo(:,:,whichOdor,i)));
  end
  [foo,Idc] = sort(dc);
  whichShufflesToPlot = Idc(1:3);
  cols = interp1([0;1],[name2rgb('darkolivegreen3'); name2rgb('darkgreen')],linspace(0,1,numShufflesToPlot));

  for i = 1:numShufflesToPlot
    Zsh = UshCbo(:,:,whichOdor,whichShufflesToPlot(i));
    coefSh  = pc'*Zsh;
    lgnd(2+i) = plot3(coefSh(1,:),coefSh(2,:),coefSh(3,:),'o-','Color',cols(i,:),'LineWidth',lineWidth,'MarkerSize',3); 
    plot3(coefSh(1,1),coefSh(2,1),coefSh(3,1),'o','Color',cols(i,:),'MarkerFaceColor',cols(i,:),'MarkerSize',5);
  end
  text(coefSh(1,1), coefSh(2,1), coefSh(3,1), 'Rec. using shuff. KCs', 'Color',cols(numShufflesToPlot,:),'FontSize',12);
  
  pcaView = load(fullfile(currDir, 'pcaView.mat'));
  view([pcaView.az pcaView.el]);
  axis tight;

  grid on;
  xlabel('PC1','FontSize',14);
  ylabel('PC2','FontSize',14);
  zlabel('PC3','FontSize',14);
end

if (any(whichPanels == 2))
  % Here we plot the distances between the PNs and the reconstructed
  % PNs, and compare them to shuffles.
  %
  data     = load(fullfile(figDir, currDir, dataDir, 'dataForPnRec.mat'));
  medSsen  = data.medSsen;
  lbSsen   = data.lbSsen;
  ubSsen   = data.ubSsen;

  medCc  = data.medCc;
  lbCc   = data.lbCc;
  ubCc   = data.ubCc;

  colors = [name2rgb('red');name2rgb('darkolivegreen3');name2rgb('gray')];  
  sfigure(FindFigureCreate('Figure 7K: PN trajectory reconstruction statistics')); clf; 
  set(gcf,'Color',[1 1 1],'Resize','off','MenuBar','none','toolbar','none', 'NumberTitle','off');
  ResizeFigure(gcf,8,8,'inches');

  Q = ComputeSubplotPositions(1,2,[],0.09,0.09,0.02,0.05,0.02,0);
  subplotp(Q,1);

  for i = 1:3
    col = colors(i,:);
    lgnd(i) = line([i i],[lbSsen(i) ubSsen(i)],'LineWidth',1,'Color',col); hold on;
    plot(i,medSsen(i),'o','Color',col,'MarkerSize',5,'MarkerFaceColor',col,'MarkerEdgeColor',col); 
  end
  xlim([0 4]);
  ylim([0 1]);
  box on;
  set(gca,'xtick',[],'xticklabel',[],'ytick',[0:0.25:1],'yticklabel',{'0','','0.5','','1'},'FontSize',12);
  title('SSE/SST', 'FontSize', 16);
  
  % Make a custom legend, because the default one is too big
  lgndY0 = 0.15;
  lgndDy = 0.05;
  lgndX0 = 0.2;
  lgndY = lgndY0;
  lgndText={'KCs, PNs','Shuff. KCs, PNs','KCs, Shuff. PNs'};
  for i = 1:3
    line([lgndX0 lgndX0+0.5],[lgndY lgndY],'Color',colors(i,:));
    text(lgndX0+0.6,lgndY-lgndDy*0.4,lgndText{i},'FontSize',6,'VerticalAlignment','Bottom');
    lgndY = lgndY - lgndDy;
  end

  % Plot the correlation coefficient data
  subplotp(Q,2);
  
  for i = 1:3
    col = colors(i,:);
    line([i i],[lbCc(i) ubCc(i)],'LineWidth',1,'Color',col); hold on;
    plot(i,medCc(i),'o','Color',col,'MarkerSize',5,'MarkerFaceColor',col,'MarkerEdgeColor',col); 
  end
  set(gca,'yaxislocation','right');
  xlim([0 4]);
  ylim([0 1]);
  box on;
  set(gca,'xtick',[],'xticklabel',[],'ytick',[0:0.25:1],'yticklabel',{'0','','0.5','','1'},'FontSize',12);
  title('Corr. Coef.', 'FontSize', 16);
end

%% Mean subtracted PCA projections
if (any(whichPanels == 3))
  % Showing the PCA reconstructions
  data = load(fullfile(figDir, currDir, dataDir, 'data100_t2_1_to_3_1.mat'));
  conn = load(fullfile(figDir, currDir, dataDir, 'prunedWeightsAdaptiveLasso.mat'));
  
  B = conn.Results(:,1:206,1)';
  B0= conn.Results(:,207,1)';
  U = data.U;
  V = data.V;
  [numSamples,numPns] = size(U);
  [numSamples,numKcs] = size(V);

  % Compute the principal components of the data
  [Upc, Ucoef] = princomp(U);
  
  % Reconstruct the trajectories
  [Urec,Upn] = ReconstructPnTrajectories(V,U,B,B0);
  
  % Mean subtract before PC projection
  sz = size(Urec);
  
  Urec1 = reshape(Urec, size(Urec,1), [])'; % timebins x cells
  Upn1  = reshape(Upn,  size(Upn, 1), [])';
  
  Urec1ms = bsxfun(@minus, Urec1, mean(Urec1,1));
  Upn1ms  = bsxfun(@minus, Upn1,  mean(Upn1, 1));
  
  Zrec     = reshape(Urec1ms', sz);
  Zpn      = reshape(Upn1ms',  sz);

  whichOdor = 25;
  Zpn  = Zpn(:,:,whichOdor);
  Zrec = Zrec(:,:,whichOdor);
  
  % Project into PC space
  coefPn  = Upc'*Zpn;
  coefRec = Upc'*Zrec;
  
  sfigure(FindFigureCreate('Figure 7J: PN trajectory reconstructions')); clf;
  set(gcf,'Color',[1 1 1],'Resize','off','NumberTitle','off'); ResizeFigure(gcf,9,7.5,'inches');

  lineWidth = 1;
  
  lgnd = [];
  lgnd(1) = plot3(coefPn(1,:),coefPn(2,:),coefPn(3,:),'bo-','LineWidth',lineWidth,'MarkerSize',3); hold on;
  plot3(coefPn(1,1),coefPn(2,1),coefPn(3,1),'o','Color','b','MarkerFaceColor','b','MarkerSize',5);
  text(coefPn(1,1), coefPn(2,1), coefPn(3,1), 'PN trajectory', 'Color','b','FontSize',12);
  
  lgnd(2) = plot3(coefRec(1,:),coefRec(2,:),coefRec(3,:),'ro-','LineWidth',lineWidth,'MarkerSize',3); 
  plot3(coefRec(1,1),coefRec(2,1),coefRec(3,1),'o','Color','r','MarkerFaceColor','r','MarkerSize',5);
  text(coefRec(1,1), coefRec(2,1), coefRec(3,1), 'Rec. using KCs', 'Color','r','FontSize',12);
  
  dataKcSh = load(fullfile(figDir, currDir, dataDir, 'lassoWeightsKcShuffle.mat'));
  [Vsh,VshCbo,Ush,UshCbo,Bsh,b0sh] = CollectDataForShuffledKcs(dataKcSh.Results);  
  numShuffles = size(Vsh,3);

  numShufflesToPlot = 3;
  % Find the shuffles which have the highest correlation for this odor.
  dc = zeros(1,numShuffles);
  for i = 1:numShuffles
    dc(i) = min(corrdist(Urec(:,:,whichOdor),UshCbo(:,:,whichOdor,i)));
  end
  [foo,Idc] = sort(dc);
  whichShufflesToPlot = Idc(1:3);
  cols = interp1([0;1],[name2rgb('darkolivegreen3'); name2rgb('darkgreen')],linspace(0,1,numShufflesToPlot));

  for i = 1:numShufflesToPlot
    Zsh       = UshCbo(:,:,:,whichShufflesToPlot(i));
    sz        = size(Zsh);
    Zsh       = reshape(Zsh, sz(1), []);
    Zshms     = bsxfun(@minus, Zsh, mean(Zsh,2));
    Zshms     = reshape(Zshms, sz);
    coefSh    = Upc'*Zshms(:,:,whichOdor);
    lgnd(2+i) = plot3(coefSh(1,:),coefSh(2,:),coefSh(3,:),'o-','Color',cols(i,:),'LineWidth',lineWidth,'MarkerSize',3); 
    plot3(coefSh(1,1),coefSh(2,1),coefSh(3,1),'o','Color',cols(i,:),'MarkerFaceColor',cols(i,:),'MarkerSize',5);
  end
  text(coefSh(1,1), coefSh(2,1), coefSh(3,1), 'Rec. using shuff. KCs', 'Color',cols(i,:),'FontSize',12);

  load(fullfile(currDir,'pcaView.mat'));
  view([az el]);
  axis tight;

  grid on;
  xlabel('PC1','FontSize',14);
  ylabel('PC2','FontSize',14);
  zlabel('PC3','FontSize',14);
end
