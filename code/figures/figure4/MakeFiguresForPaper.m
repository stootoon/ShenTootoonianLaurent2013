function MakeFiguresForPaper(whichPanels, varargin)
% function MakeFiguresForPaper(whichPanels, ['dataDir'='originalData'])
%
% whichPanels = 1 => Figure 4A: Example fit procedure
% whichPanels = 2 => Figure 4B/S4A-C: Model fit grids, r2, snr
% whichPanels = 3 => Figure S4D: The model distribution
% whichPanels = 4 => Figure S4E: Scaling factor vs. mixture level
% whichPanels = 5 => Figure S4F: SNR Angle results
% whichPanels = 5 => Figure S4G: Weight ratio

p = inputParser;
p.addOptional('dataDir', 'originalData');
p.parse(varargin{:});

whichFigure = 4;
figDir  = GetDataDirForFigure(whichFigure);
currDir = ''; % We're not in a subdirectory relative to the figDir
dataDir = p.Results.dataDir;

pnDataFile = 'pnComplexMixtures0_0to6_0_50ms.mat';
allCmps    = 'ABCDWXYZ';
MUSDSE     = @(x) [mean(x) std(x) std(x)/sqrt(numel(x))];

if (any(whichPanels==1))
  % 4A: Example fit procedure

  % Load the data
  resultsFile = 'bestModelsForMixturesResults_allOdors_sigmaReg_1_0_sigmaLag_1_0_lagLimit_3_vAlpha_1_0_vBeta_1_0.mat';
  RS = load(fullfile(figDir, currDir, dataDir, resultsFile));
  fitOptions = {'sigmaLag',RS.sigmaLag,'sigmaReg',RS.sigmaReg,'lagLimit',RS.lagLimit,'variancePriorAlpha',RS.variancePriorAlpha, 'variancePriorBeta', RS.variancePriorBeta, 'whichBinsToFit', RS.whichBinsToFit};
  
  % Build the data for the reconstruction
  whichCell = 6;
  whichOdor = 'odorABCD';
  cmpsInOdor = whichOdor(5:end);
  numCmpsInOdor = numel(whichOdor) - 4;
  whichOdorInd  = odor_name_to_index(whichOdor);
  
  % Do free models adding 1 component at a time, then do the best model.
  WhichModelAndInputConfig = [];
  WhichModelAndInputConfig(1).model  = 1;
  WhichModelAndInputConfig(1).config = 0*ones(1, numCmpsInOdor);
  
  for i = 1:numCmpsInOdor
    WhichModelAndInputConfig(i+1).model  = 6; % Free model
    WhichModelAndInputConfig(i+1).config = [ones(1,i) zeros(1, numCmpsInOdor-i)]; % Free model
  end
  WhichModelAndInputConfig(i+2).model  = RS.Results{whichCell, whichOdorInd-12}{1};
  WhichModelAndInputConfig(i+2).config = RS.Results{whichCell, whichOdorInd-12}{2};

  indCmps = GetComponentsForOdor(whichOdorInd);
  Data = load(fullfile(figDir, currDir, dataDir, pnDataFile));
  Xobs = squeeze(Data.pnCbot(whichCell,:,indCmps+3,:));
  yobs = squeeze(Data.pnCbot(whichCell,:,whichOdorInd,:));
  
  Ye = [];
  descr = {};
  % Fit all of the models
  for i = 1:numel(WhichModelAndInputConfig)
    whichModel       = WhichModelAndInputConfig(i).model;
    whichInputConfig = WhichModelAndInputConfig(i).config;
    F = FitAllModelsForObservations2LaplaceGeneral(Xobs, yobs, fitOptions{:}, 'fitMode', 'single', 'whichInputConfig', whichInputConfig, 'whichModel', whichModel);
    Ye(i,:) = F.Model.test(F.Model.X);
    str = F.descr;
    fprintf('Model: %s\nR2: %1.2e\n\n', str, F.Model.r2);
    for j = 1:numCmpsInOdor
      str = strrep(str,['x_' num2str(j)], cmpsInOdor(j));
    end
    descr{i} = str;
  end

  % Plot the different fits
  sfigure(FindFigureCreate('Figure 4A: Example fit procedure')); clf;
  set(gcf,'Color',[1 1 1],'Resize','off','NumberTitle','off');
  ResizeFigure(gcf,15,6, 'inches');
  Q = ComputeSubplotPositions(5,3,{...
      {1,1}, {2,1}, {3,1}, {4,1}, {5,1},...
      {1,2}, {1,3},...
      {3,2}, {3,3},...
      {5,2}, {5,3},...
                   },...
                              0.025,0.01,0.025,0.075,0.1,0.01);
  % Plot the observation
  xl = [Data.tfit(1) Data.tfit(end)];
  
  yl = [-0.1 3];
  ylCmp = [-0.1 1.5];
  
  ytick  = 0:0.5:3;
  
  patchy = [-0.09 2.9 2.9 -0.09];
  patchyCmp = [-0.09 1.45 1.45 -0.09];
  
  patchCol = name2rgb('gray90');

  subplotp(Q,1);
  patch([2 2 2.5 2.5], patchy, patchCol, 'EdgeColor', 'none'); hold on;   
  h1 = plot(Data.tall, mean(yobs,2),'k','LineWidth',1);
  set(legend(h1, [cmpsInOdor '(t)']),'box','off','color','none','FontSize',10);
  xlim(xl);
  ylim(yl);
  set(gca,'xticklabel',[],'ytick', ytick, 'yticklabel',[]);
  box on;
  % Plot the components
  cmpCols = GetComponentColors;
  for i = 1:numCmpsInOdor
    subplotp(Q,i+1);
    patch([2 2 2.5 2.5], patchyCmp, patchCol, 'EdgeColor', 'none'); hold on;   
    h1 = plot(Data.tall, mean(Xobs(:,i,:),3),'Color',cmpCols(indCmps(i),:),'LineWidth',1);
    set(legend(h1,[allCmps(indCmps(i)) '(t)']),'box','off','color','none','FontSize',10);    
    ylim(ylCmp);
    xlim(xl);
    set(gca,'xticklabel',[],'yticklabel',[]);
    box on;
  end
  
  % Now plot the fits
  fitColors = [name2rgb('gray70');
               name2rgb('skyblue2');
               name2rgb('dodgerblue');
               name2rgb('dodgerblue2');
               name2rgb('dodgerblue3');
               name2rgb('magenta');];
  for i = 1:size(Ye,1)
    subplotp(Q,5+i);
    patch([2 2 2.5 2.5], patchy, patchCol, 'EdgeColor', 'none'); hold on;   
    plot(Data.tall, mean(yobs,2), 'k', 'LineWidth',1);    
    plot(Data.tfit, Ye(i,:),      'Color', fitColors(i,:), 'LineWidth',1);
    ylim(yl);
    xlim(xl);
    set(gca,'xticklabel',[],'yticklabel',[],'ytick',ytick);
    box on;
  end

  % Now annotate the plots
  subplotp(Q,5);
  set(gca,'xtick',[0:3]+2,'xticklabel',arrayfun(@num2str, [0:3],'UniformOutput',false));
  set(gca,'ytick',[0:0.5:1.5], 'yticklabel',{'0','','1',''});
  xlabel('Time (s)', 'FontSize', 12);

  subplotp(Q,1);
  yt = [0:0.5:3];
  ytl= arrayfun(@(x) '', yt, 'UniformOutput',false);
  ytl{1}   = '0';
  ytl{end} = '3';
  set(gca,'ytick', yt, 'ytickLabel', ytl);
  h  = title('DATA','FontSize',16,'FontWeight','bold', 'HorizontalAlign','Left');
  pos = get(h,'Position');
  pos(1) = min(xlim);
  pos(2) = 1.2*max(ylim);
  set(h,'Position',pos);

  subplotp(Q,3);
  ylabel('FR/20 (Hz)', 'FontSize', 12);

  subplotp(Q,6);
  h  = title('FITS','FontSize',16,'FontWeight','bold', 'HorizontalAlign','Left');
  bigTitle = copyobj(h, get(h,'Parent'));

  pos = get(bigTitle,'Position');
  pos(1) = min(xlim);
  pos(2) = 1.2*max(ylim);
  set(bigTitle,'Position',pos);

  title('Constant Model','FontSize', 14, 'HorizontalAlign', 'Center');
  
  subplotp(Q,7);
  title('Only A','FontSize', 14, 'HorizontalAlign', 'Center');
  
  subplotp(Q,8);
  title('A and B','FontSize', 14, 'HorizontalAlign', 'Center');
  
  subplotp(Q,9);
  title('A and B and C','FontSize', 14, 'HorizontalAlign', 'Center');

  subplotp(Q,10);
  title('A and B and C and D','FontSize', 14, 'HorizontalAlign', 'Center');
  
  subplotp(Q,11);
  title('Unit-scaled A and B and C*','FontSize', 14, 'FontWeight','bold', 'HorizontalAlign', 'Center');

  for i = 1:6
    subplotp(Q,i+5);
    xlabel(strrep(descr{i}, 'y(t)', 'fit(t)'));
  end
end

%% 4B/SA-C: Model fits, r2, snr
if (any(whichPanels==2))
  % 4B/SA-C: Model fits, r2, snr

  dataForFigure = load(fullfile(figDir, currDir, dataDir, 'dataForFigure'));
  resultsFile   = 'bestModelsForMixturesResults_allOdors_sigmaReg_1_0_sigmaLag_1_0_lagLimit_3_vAlpha_1_0_vBeta_1_0.mat';
  ResultsStruct = load(fullfile(figDir, currDir, dataDir, resultsFile));
  
  snr = ComputeSnrForAllCellsAndMixtures('dataFile',fullfile(figDir, currDir, dataDir, 'pnComplexMixtures0_0to13_0_50ms.mat')); % ComputeSnrForAllCellsAndMixtures();
  snr = snr(:,1:size(ResultsStruct.Results,2));

  snr(isinf(snr)) = max(snr(~isinf(snr)));
  snrDb           = 10*log10(snr); % whereever snr is zero it will be -Inf.
  snrDb(isinf(snrDb)) = min(snrDb(~isinf(snrDb(:)))); 

  numCells    = 174;
  cmFits      = reshape(dataForFigure.colorMaps{1},numCells,[],3);
  numMixtures = size(cmFits,2);

  M = reshape(1:numCells*numMixtures, numCells, numMixtures);
  
  aspectRatio = numMixtures/numCells;

  pnGroupSizes = cellfun(@numel, dataForFigure.I1);
  pnGroup1     = vertcat(dataForFigure.I1{:});
  pnGroupMix   = [vertcat(dataForFigure.Imix{:});vertcat(dataForFigure.Iamb{:})];
  
  sfigure(FindFigureCreate('Figure 4B/S4A-C: Single Pn Mixture Response Fits')); clf;
  set(gcf,'Color',[1 1 1],'Resize','off','Units','normalized','NumberTitle','off');
  ResizeFigure(gcf,6,12,'inches');

  nrLgnd     = 2;
  nrLgndGap  = 2;
  
  nrMixGap   = 10;
  nrMixLgnd  = 8;
  
  nr = 8*(nrLgnd+2*nrLgndGap) + numel(pnGroup1) + nrMixGap + nrMixLgnd + numel(pnGroupMix);
  nc = 3;

  % Generate the subplot positions for the response type plots.  The
  % layouts for these is non-standard. The codes below indicate the
  % layout for each column of the plot. The codes are:
  %
  % l#: a legend for the corresponding component (1 = A, 2 = B, etc...) 
  % #: Grid for the cells with #-type responses 
  % g: gap
  % gmix: bigger gap to separate the single component results from the mixture.
  % lmix: Legend for the mixture/ambient-type response.
  % _l: Gap the size of a legend.
  % _lmix: Gap the size of a mixture legend.
 
  columnLayout = {{'l1','g',1,'g','l2','g',2,'g','l3','g',3,'g','l4','g',4,'g',...
                   'l5','g',5,'g','l6','g',6,'g','l7','g',7,'g','l8','g',8,'gmix',...
                   'lmix','g',9},...
                  {'l1','g',1,'g','_l','g',2,'g','_l','g',3,'g','_l','g',4,'g',...
                   '_l','g',5,'g','_l','g',6,'g','_l','g',7,'g','_l','g',8,'gmix',...
                   '_lmix','g',9},...
                  {'l1','g',1,'g','_l','g',2,'g','_l','g',3,'g','_l','g',4,'g',...
                   '_l','g',5,'g','_l','g',6,'g','_l','g',7,'g','_l','g',8,'gmix',...
                   '_lmix','g',9}};

  plotPositions = {};
  ip = 1;
  lgnds  = arrayfun(@(i) ['l' num2str(i)], 1:8, 'UniformOutput', false);
  
  % Now load plotPositions according to the layout above.
  for i = 1:3
    currRow = 1;
    layout = columnLayout{i};
    for j = 1:numel(layout)
      switch layout{j}
       case lgnds
        plotPositions{ip} = {currRow+(0:nrLgnd),i};
        currRow = currRow + nrLgnd;
        ip = ip+1;
       case 'lmix'
        plotPositions{ip} = {currRow+(0:nrMixLgnd),i};
        currRow = currRow + nrMixLgnd;
        ip = ip+1;        
       case '_l'
        currRow = currRow + nrLgnd;
       case '_lmix'
        currRow = currRow + nrMixLgnd;
       case 'g'
        currRow = currRow + nrLgndGap;
       case 'gmix'
        currRow = currRow + nrMixGap;
       case {1,2,3,4,5,6,7,8}
        numCells = pnGroupSizes(layout{j});
        plotPositions{ip} = {currRow+(0:numCells-1), i};
        ip = ip+1;
        currRow = currRow + numCells;
       case 9
        numMix  = numel(pnGroupMix);
        plotPositions{ip} = {currRow + (0:numMix-1), i};
        ip = ip+1;
        currRow = currRow + numMix;
       otherwise
        error('Unknown column signature "%s".', layout{j});
      end
    end
  end

  % Setup the plots
  Q = ComputeSubplotPositions(nr,nc,plotPositions, 0.075,0.025,0.05,0.05,0.01,0.0);

  B = GetOdorNamesAsBinaryVectors('full');
  B = B(first(13:size(B,1),numMixtures),:);
  
  Bcol = B*diag(1:8);
  
  mixtureLevel = sum(B,2);
  
  cmpCols  = GetComponentColors;
  whichCms = [1 4 3];
  currPlot = 1;
  cmpNames = 'ABCDWXYZ';
  for cm = 1:3
    rowOffset = 0;
    indAll = 1:size(M,2);
    layout = columnLayout{cm};

    if (cm==2)
        colLegend = interp1([0; 0.1; 0.25; 0.5; 0.75; 0.9; 1.0],...
                            [name2rgb('black');...
                            name2rgb('red');...
                            name2rgb('orange');...
                            name2rgb('white');...
                            name2rgb('deepskyblue2');...
                            name2rgb('royalblue3'); 
                            name2rgb('royalblue3')],linspace(0,1,128));      
    elseif (cm==3)
        colLegend = interp1([min(snrDb(:)); -5; 0; 3; 5; 10; 15; 20; max(snrDb(:))],...
                            [name2rgb('royalblue');...
                            name2rgb('royalblue');...
                            name2rgb('black');...                
                            name2rgb('darkred');...
                            name2rgb('red');...
                            name2rgb('yellow');...
                            name2rgb('white');...
                            name2rgb('white');...
                            name2rgb('white');],linspace(-5,20,128));
    end
    
    % Now plot the actual grids
    for i = 1:numel(layout)
      switch layout{i}
       case lgnds % Plot the legends
        if (cm==1)
          whichCmp = str2num(layout{i}(2));
          indCmp   = intersect(indAll, find(B(:,whichCmp))); 
          mlCmp    = mixtureLevel(indCmp);
          indNoCmp = setdiff(indAll, indCmp);
          mlNoCmp  = mixtureLevel(indNoCmp);
          colLegend= [ones(size(indCmp(:)));zeros(size(indNoCmp(:)))]*cmpCols(whichCmp,:);
        end
        subplotp(Q, currPlot);
        imagesc(1:size(colLegend,1),[1 size(colLegend,1)]);
        set(gca,'xtick',[],'ytick',[]); 
        if (isequal(layout{i}, 'l1'))
          switch cm
           case 1
            xlabel('Odors','FontSize',14);
            set(gca,'XAxisLocation','top');
           case 2
            xl = [0 1];
            imagesc(xl, [1 nrLgnd], linspace(xl(1), xl(2),100));            
            xlabel('R^2','FontSize',14);
            set(gca,'xtick',[0 1],'ytick',[]);
            set(gca,'XAxisLocation','top');
           case 3
            xl = [-5 20];
            imagesc(xl, [1 nrLgnd], linspace(xl(1), xl(2),100));
            set(gca,'xtick',[0 3 20],'ytick',[]);
            xlabel('SNR (dB)','FontSize',14);
            set(gca,'XAxisLocation','top');
          end            
        end

        colormap(colLegend);
        freezeColors;
        currPlot = currPlot + 1;
       case {1,2,3,4,5,6,7,8} % Plot the single-component dominated responses
        whichCmp = layout{i};
        indCmp   = intersect(indAll, find(B(:,whichCmp))); 
        mlCmp    = mixtureLevel(indCmp);
        indNoCmp = setdiff(indAll, indCmp);
        mlNoCmp  = mixtureLevel(indNoCmp);
   
        subplotp(Q, currPlot);
        imagesc(M(sum(pnGroupSizes(1:whichCmp-1))+(1:pnGroupSizes(whichCmp)),[indCmp(:);indNoCmp(:)]),[1 numel(M)]);
        set(gca,'xtick',[],'ytick',[]); 
        if (cm==1)
          ylabel(cmpNames(whichCmp),'FontSize',12);
        end
        axis equal; axis tight;
        colormap(dataForFigure.colorMaps{whichCms(cm)});    
        freezeColors;
        currPlot = currPlot + 1;
       case 9 % Plot the mixture/ambient type responses
        subplotp(Q, currPlot);
        imagesc(M(sum(pnGroupSizes)+(1:numel(pnGroupMix)),:),[1 numel(M)]);    
        colormap(dataForFigure.colorMaps{whichCms(cm)});
        set(gca,'xtick',[],'ytick',[]);
        if (cm==1)
          ylabel('Mix/Amb','FontSize',12);
        end
        axis equal; axis tight;
        freezeColors;
        currPlot = currPlot + 1;
       case 'lmix' % Plot the mixture legend
        subplotp(Q, currPlot);
        imagesc(Bcol',[0 8]);
        colormap([0 0 0; cmpCols]);
        axis equal; axis tight;
        set(gca,'xtick',[],'ytick',[]); 
        freezeColors;
        currPlot = currPlot + 1;
       case 'g' % All the rest are gaps
       case 'gmix'
       case '_l'
       case '_lmix'
       otherwise
        error('Unknown layout element "%s".', layout{i});
      end
    end    
  end
end

if (any(whichPanels==3))
  % Make Figure S4D, the Model distributions.
  ModelStats = LoadVarFromMatFileByName(fullfile(figDir, currDir, dataDir, 'modelStats'),'ModelStats');
    
  modelTypes = zeros(8,4); % 8 components, 4 types: unit, unit-lagged, scaled, scaled-lagged
  ilag       = [ModelStats.ilag];
  itype      = [ModelStats.itype];
  istruct    = [ModelStats.istruct];
  
  indSingle = find(itype>=1 & itype<9);
  subs = [itype(indSingle)' istruct(indSingle)'+2*(ilag(indSingle)>0)'];
  singleModelDistribution = accumarray(subs,1);
  singleModelDistribution = singleModelDistribution(:,[1 3 2 4]); % 1,1-lag, k,k-lag
  
  indMix = find(itype == 10);
  subs   = [istruct(indMix)'+3*(ilag(indMix)>0)'];
  mixModelDistribution = accumarray(subs,1)';
  mixModelDistribution = mixModelDistribution([1 4 2 5 3 6]); %1,1-lag, k,k-lag, free,free-lag
  
  % Make the stacked chart
  sfigure(FindFigureCreate('Figure S4D: Model distribution'));
  clf; set(gcf,'Color',[1 1 1],'Resize', 'off', 'NumberTitle','off');
  ResizeFigure(gcf,8,3,'inches');
  data = zeros(9,6);
  data(1:8,1:4) = bsxfun(@rdivide, singleModelDistribution, sum(singleModelDistribution,2));
  data(9,   :)  = mixModelDistribution/sum(mixModelDistribution);
  bar([1:8 9.5], data, 0.6, 'stacked');
  cm = [name2rgb('dodgerblue3');
        name2rgb('dodgerblue1');
        name2rgb('sgiteal');
        name2rgb('sgichartreuse');
        name2rgb('red');
        name2rgb('darkorange')];  
  colormap(cm);
  xlim([0.5 14]);
  ylim([0 1]);
  cmps = 'ABCDWXYZ';
  set(gca,'xtick',[1:8 9.5], 'xticklabel', [arrayfun(@(x) {x}, cmps) 'Mix'], 'ytick',[0:0.2:1],'ticklength',[0 0],'FontSize',12);
  xlabel('Response Type','FontSize',       14);
  ylabel('Frequency Observed', 'FontSize', 14);
  legend('Unit','Unit (lag)', 'Scaled', 'Scaled (lag)', 'Free', 'Free (lag)','Location','East');
  TightenAxesToFigure;
end

if (any(whichPanels==4))
  ModelStats    = LoadVarFromMatFileByName(fullfile(figDir, currDir, dataDir, 'modelStats'),'ModelStats');
  dataForFigure = load(fullfile(figDir, currDir, dataDir, 'dataForFigure'));
  B  = GetOdorNamesAsBinaryVectors('full');
  B  = B(13:end,:);
  ml = sum(B,2);

  istruct = [ModelStats.istruct];
  istruct = reshape(istruct, size(ModelStats));
  istructIsScaled = (istruct == 1 | istruct==2); % Look for the scaled models
  counts = zeros(8,4); % 8 components, 4 differences (3-2,4-2,5-2,8-2)
  diffs  = zeros(8,4);
  mls  = [2:5 8];
  recs = [];

  sfigure(FindFigureCreate('Figure S4E: Scaling factor vs. mixture level'));
  clf; set(gcf,'Resize','off','Color',[1 1 1],'NumberTitle','off');
  ResizeFigure(gcf,8,3,'inches');
  Q = ComputeSubplotPositions(2,4,[],0.075,0.01,0.01,0.01,0.2,0.02);
  cmpColors = GetComponentColors;
  cmps = 'ABCDWXYZ';
  for i = 1:8 % Loop over components
    subplotp(Q,i);
    inds = dataForFigure.I1{i};
    for j = 1:numel(inds) % Loop over cells that produced responses of this type
      newRec = nan(1,5);
      for iml = 1:5    % Find the mean k at each ml
        thisMl = mls(iml);
        indMl = find(ml==thisMl & istructIsScaled(inds(j),:)'); % Get all the responses that were scaled or unit.
        if (~isempty(indMl)) % Found some at this ml
          w1 = zeros(1,numel(indMl));
          for k = 1:numel(indMl) % Grab the weights
            w = ModelStats(inds(j),indMl(k)).weights;
            w = unique(w(~isnan(w)));
            assert(~isempty(w) && numel(w)==1,'Expected exactly 1 weight for scaled model.');
            w1(k) = w;
          end
          newRec(iml) = mean(w1); % Store their mean value
        end
      end

      cols = cmpColors(i,:)*(rand*0.4+0.6);
      indNonNan = find(~isnan(newRec));    
      if (numel(indNonNan)>2)
        plot(mls(indNonNan),newRec(indNonNan), 'o-','Color',cols,'LineWidth',0.5,'MarkerFaceColor',cols,'MarkerSize',3); hold on;
      end
      recs(size(recs,1)+1,:) = [i j newRec];
    end
    xlim([1.5 8.5]);
    ylim([0 1.5]);
    set(gca,'xtick',[2 3 4 5 8],'ytick',[0:0.5:1.5],'FontSize',12);
    if (i==5)
      % W; bottom left, will hold the annotatations
      xlabel('Mixture Level',  'FontSize', 14);
      ylabel('Scaling Factor', 'FontSize', 14);
    else
      set(gca,'xticklabel',[],'yticklabel',[]);        
    end
    
    text(8, 1.4, [cmps(i) '-type'],'VerticalAlignment','top','HorizontalAlignment','right','FontSize',12);
  end
end

if (any(whichPanels==5))
%% Plot the SNR angle results
  snrAngleData = load(fullfile(figDir, currDir, dataDir, 'snrAngleData'));
  thAll  = snrAngleData.thAll;
  thPref = snrAngleData.thPref;
  
  sfigure(FindFigureCreate('S4F: SNR angle'));
  clf; set(gcf,'Color',[1 1 1],'Resize','off','NumberTitle','off');
  ResizeFigure(gcf,8,6,'inches');
  bins = -pi:pi/128:pi;
  bw = 0.12/2;
  [f1,x1,bw1] = ksdensity(thAll,bins,'width',bw);
  [f2,x2,bw2] = ksdensity(thPref,bins,'width',bw);
  binSize = first(diff(bins));
  
  h1 = area(x1,f1*binSize); hold on;
  set(h1,'FaceColor',name2rgb('royalblue'),'EdgeColor', name2rgb('royalblue')); 
  set(get(h1,'Children'),'FaceAlpha', 0.2);

  h2 = area(x2,f2*binSize); 
  set(h2,'FaceColor',name2rgb('red'),      'EdgeColor', name2rgb('red')); 
  set(get(h2,'Children'),'FaceAlpha', 0.2);
  
  legend([h1 h2], 'All 1-Comp. Fits', 'Preferred-Comp. Fits Only','Location','NorthEast');
  
  set(gca,'xtick',-pi/2:pi/2:pi,'xticklabel',{'-pi/2','0','pi/2','pi'},'ytick',[0:0.01:0.03],'FontSize',12,'tickLength',[0.005 0.005]);
  xlabel('SNR Angle', 'FontSize', 14);
  ylabel('Frequency Observed','FontSize', 14);
  TightenAxesToFigure;
  xlim([-pi/2 pi]);
end

if (any(whichPanels==6))
  cutoff = 0.2;
  data   = load(fullfile(figDir, currDir, dataDir, 'weightRatioData'));
  weightRatio = data.weightRatio;
  
  sfigure(FindFigureCreate('Figure S4G: Weight ratio')); clf;
  set(gcf,'Color',[1 1 1],'Resize','off','NumberTitle','off');
  ResizeFigure(gcf, 8, 4,'inches');
  
  [f,xi] = ksdensity(weightRatio,linspace(-2,4,100),'width',1/20); % There isn't much beyond 4.
  indBefore = find(xi<=cutoff);
  indAfter  = find(xi>=cutoff);
  xib = xi(indBefore);
  yib = f(indBefore)*first(diff(xi));
  hb = area(xib,yib); 
  set(hb,'FaceColor',name2rgb('orangered'),'EdgeColor',name2rgb('orangered'),'LineWidth', 2); hold on;
  set(get(hb,'Children'),'FaceAlpha',0.2);
  
  xia = xi(indAfter);
  yia = f(indAfter)*first(diff(xi));
  ha  = area([xib(end) xia], [yib(end) yia]);
  set(ha,'FaceColor',name2rgb('goldenrod'),'EdgeColor',name2rgb('goldenrod'),'LineWidth',2);
  set(get(ha,'Children'),'FaceAlpha',0.2);

  xlim([-2 3]);
  ylim([0 0.08]);
  VerticalLine(gca,cutoff,'Color', name2rgb('gray50'),'LineStyle', '--');
  set(gca,'xtick',[-2:3],'ytick',0:0.02:0.08,'FontSize',12);
  xlabel('Weight Ratio','FontSize',14);
  ylabel('Frequency Observed', 'FontSize', 14);
  text(-1,0.07,'Suppression','Color',name2rgb('orangered'),'VerticalAlignment','top','HorizontalAlignment','center','FontSize',14);
  text( 1,0.05,'Redundancy', 'Color',name2rgb('goldenrod'),'VerticalAlignment','top','HorizontalAlignment','center','FontSize',14);
  TightenAxesToFigure;
end