function MakeFiguresForPaper(whichFigures, varargin)
% MakeFiguresForPaper(whichFigures, varargin)
%
% Helper function that makes most of the panels for Figure 2 and
% Figure S2. The figure IDs are:
%
% 1: Preferred response types, Figure 2F, S2A.
% 2: R^2 values, Figure S2B.
% 3: SNR values, Figure S2C.
% 4: Example fit procedure, Figure 2E.
% 5: Model distribution pie charts, Figure S2D
% 6: Effect of mixture on fits, Figure S2E.
% 7: SNR angles for each response type (Figure S2F).

p = inputParser;
p.addOptional('labelFontSize', 8);
p.addOptional('axisFontSize', 6);
p.addOptional('doExport', 0);
p.addOptional('dataDir', 'originalData');
p.parse(varargin{:});

addpath(fileparts(mfilename('fullpath')));

labelFontSize = p.Results.labelFontSize;
axisFontSize  = p.Results.axisFontSize;

dataDir       = p.Results.dataDir;
figDir        = GetDataDirForFigure(2);
thisDir       = GetCurrentDirFromPathString(fileparts(mfilename('fullpath')));

XLABEL_ = @(x) xlabel(x,'FontSize',labelFontSize);
YLABEL_ = @(x) ylabel(x,'FontSize',labelFontSize);

%% The first three figures
%  1: Makes the figure showing the odor each PN likes
%  2: Makes the figure showing the R2's for each non-constant response
%  3: Makes the figure showing the SNRs for each non-constant response

for i = 1:numel(whichFigures)
  if (whichFigures(i)>=1 && whichFigures(i)<=3)

    dataForFigure = load(fullfile(figDir, thisDir, dataDir, 'dataForFigure'));
    pnGroupSizes = cellfun(@length, dataForFigure.I);
    
    cm = dataForFigure.colorMaps{whichFigures(i)+3}; % The color map for the data
    
    [mixtureInds, mixtureVals] = GetBinaryMixtureIndsForResponseLinearity;
  
    % We'll plot Octanol since it has the most, and then citral ontop of the mixtures
    numCitral = pnGroupSizes(2);
    numOctanol= pnGroupSizes(3);
    numMixAmb = sum(pnGroupSizes(4:5));
    numLegendCols = 3;
    numSpacer     = 2;
    
    numCols = numLegendCols + numSpacer + numCitral + numSpacer + numOctanol + numSpacer + numMixAmb;
    numBinaryMixtures = numel(mixtureInds);
    numRows = 11; % We're only going to use the morphs

    switch whichFigures(i)
     case 1
      figureName = 'Figure 2F/S2A: Best fits';
     case 2
      figureName = 'Figure S2B: Fit quality';
     case 3
      figureName = 'Figure S2C: SNR';
    end
    
    sfigure(FindFigureCreate(figureName)); clf;
    set(gcf,'Color',[1 1 1],'Resize','off','Units','normalized','NumberTitle','off');
    ResizeFigure(gcf,12,2,'inches');
    Q = ComputeSubplotPositions(numRows, numCols,...
                       {{1:numRows, 1:numLegendCols},...
                        {1:numRows, numLegendCols + numSpacer + (1:numCitral)},...
                        {1:numRows, numLegendCols + numCitral + 2*numSpacer +  (1:numOctanol)},...
                        {1:numRows, numLegendCols + numCitral + numOctanol + 3*numSpacer + (1:numMixAmb)}},...
                                0.05,0.01,0.0,0.2,0.01,0.0);

    M = reshape(1:sum(pnGroupSizes)*numBinaryMixtures, [], numBinaryMixtures)';
    indCitToOct = 11:-1:1; % Ordered to have citral on the bottom 
                           
    % Plot citral 
    subplotp(Q,2);
    Mcit = M(indCitToOct, sum(pnGroupSizes(1))+(1:pnGroupSizes(2)));
    imagesc(Mcit,[1 sum(pnGroupSizes)*numBinaryMixtures]);
    colormap(cm);  
    set(gca,'ytick',[],'xtick',[]);
    axis equal; axis tight;
    freezeColors;
    title('Citral PNs', 'FontSize', 14);

    % Now octanol
    subplotp(Q,3);
    Moct = M(indCitToOct, sum(pnGroupSizes(1:2))+(1:pnGroupSizes(3))); 
    imagesc(Moct,[1 sum(pnGroupSizes)*numBinaryMixtures]);
    colormap(cm);  
    set(gca,'ytick',[],'xtick',[]);
    axis equal; axis tight;
    freezeColors;
    title('Octanol PNs', 'FontSize', 14);
  
    % Now mixtures
    subplotp(Q,4);
    Mmix = M(indCitToOct, sum(pnGroupSizes(1:3))+(1:sum(pnGroupSizes(4:5))));
    imagesc(Mmix,[1 sum(pnGroupSizes)*numBinaryMixtures]);
    colormap(cm);  
    set(gca,'ytick',[],'xtick',[]);
    axis equal; axis tight;
    freezeColors;
    title({'Mix./','Ambig.','PNs'}, 'FontSize', 14);
    
    % Now the legends
    switch whichFigures(i)
     case 1
      subplotp(Q, 1);
      cmLegend = mixtureVals(1:11,:)/140*[1 0 0; 0 1 0];
      imagesc((11:-1:1)'*ones(1,numLegendCols));
      colormap(cmLegend);
      set(gca,'xtick',[],'ytick',[]);
      axis equal; axis tight;
      freezeColors;  
      ylabel('Mixture', 'FontSize', 12);
     case 2
      subplotp(Q,1);
      r2l = linspace(0,1,100);
      cmLegend = interp1([0; 0.1; 0.25; 0.5; 0.75; 0.9; 1.0],...
                         [name2rgb('black');...
                          name2rgb('red');...
                          name2rgb('orange');...
                          name2rgb('white');...
                          name2rgb('deepskyblue2');...
                          name2rgb('royalblue3'); 
                          name2rgb('royalblue3')],r2l);
      numThisCols = round(numel(r2l)/numRows*numLegendCols);
      
      imagesc([1 numThisCols], [1 numel(r2l)], r2l(:)*ones(1,numThisCols), [0 1]);
      axis xy;
      colormap(cmLegend);
      set(gca,'xtick',[],'ytick',[1 numel(r2l)],'yticklabel',{'0','1'});
      axis equal; axis tight;
      freezeColors;
      ylabel('R^2', 'FontSize', 12);
     case 3
      subplotp(Q,1);
      snrDbl = linspace(-5,20,100);
      cmLegend = interp1([-100;-10; -5; 0; 5; 10; 15; 20; 25; 100],...
                         [name2rgb('royalblue');...
                          name2rgb('royalblue');...
                          name2rgb('royalblue');...
                          name2rgb('black');...
                          name2rgb('red');...
                          name2rgb('orange');...
                          name2rgb('yellow');...
                          name2rgb('white');...
                          name2rgb('white');...
                          name2rgb('white');],snrDbl(:));
      
      numThisCols = round((max(snrDbl(:)) - min(snrDbl(:)))/numRows*numLegendCols);
      imagesc([1 numThisCols], [-5 20], snrDbl(:)*ones(1,numThisCols),[-5 20]);
      axis xy;
      colormap(cmLegend);
      set(gca,'xtick',[],'ytick',[0 3 20],'yticklabel',{'0','3','20'});
      axis equal; axis tight;
      freezeColors;
      ylabel('SNR(dB)');
      
    end
  end
end

%% Figure 2E: Example fit procedure
if (any(whichFigures==4))
  resultsFile = fullfile(figDir, thisDir, dataDir, 'resultsBestModel3Bayes1Laplace_slag100_sreg100_maxLag3_funcModelPrior.mat');
  Output = load(resultsFile);
  
  pnDataFile = fullfile(figDir, thisDir, dataDir, 'delayRegressData.mat');
  
  whichCell = 128;
  whichMixture = [120 140];
  
  mixtureInd = GetResponseLinearityIndexForBinaryMixtureConcentrationPair(whichMixture(1), whichMixture(2));    
  Fit = FitAllModelsForCellAndMixtureBayes1Laplace(whichCell, mixtureInd, 'pnCbotDataFile', pnDataFile, Output.fitOptions{:});    

  sfigure(FindFigureCreate('Figure 2E: Example fit procedure')); clf; 
  set(gcf,'Color',[1 1 1], 'Resize','off','NumberTitle','off');
  ResizeFigure(gcf, 12, 4, 'inches');

  Q = ComputeSubplotPositions(3,6,{...
      {1,1},{2,1},{3,1},...
      {1,2},{2,2},{3,2},...
      {1,3},{2,3},{3,3},...
      {1,4},{2,4},{3,4},...
      {1,5},{2,5},{3,5},...
      {1,6},{2,6},{3,6},...
                   }, 0.05,0.02,0.04,0.1,0.2,0.02);
  yl = [-0.1 1.6];
  patchy = [-0.08 1.55 1.55 -0.08];
  patchColor = name2rgb('gray90');
  octData = Fit.X(:,1);
  citData = Fit.X(:,2);
  mixData = Fit.X(:,end);
  tfit    = Fit.tfit;
  yfit    = nan(6,1)*tfit(:)';
  whichModels = [1 12 13 14 15];
  plotData = [1 1;
              0 0;
              1 0;
              0 1;
              1 1;
              1 1];
  bestLl = Fit.Models{Fit.bestModel}.logModelLikelihood;

  % Load the fits, and save the descriptions
  modelDescr = {};
  for i = 1:numel(whichModels)
    M = Fit.Models{whichModels(i)};
    A = Fit.Args{whichModels(i)};
    yfit(i+1,:) = M.test(A);
    disp(M.descr)
    modelDescr{i} = M.descr;
    
    fprintf('Relative LL: %1.3f\n',       M.logModelLikelihood - bestLl);
    fprintf('SSE: %1.3f, R2: %1.3f\n\n',  M.sse, M.r2)
  end
  
  % Now make the plots
  fitColors = [name2rgb('black');
               name2rgb('gray70');
               name2rgb('darkturquoise');
               name2rgb('magenta');
               name2rgb('cobalt');
               name2rgb('blue')];

  for i = 1:size(yfit,1)
    subplotp(Q,(i-1)*3+1);
    patch([2 2 2.3 2.3], patchy, patchColor,'EdgeColor','none'); hold on;
    plot(tfit, mixData, 'Color', 'k', 'LineWidth', 1); 
    if (all(~isnan(yfit(i,:))))
      plot(tfit, yfit(i,:),'Color', fitColors(i,:),'LineWidth',1);
    end
    xlim([tfit(1) tfit(end)]);
    ylim(yl);
    box on;
    set(gca,'xticklabel',[],'yticklabel',[]);
    
    if (plotData(i,1))
      subplotp(Q,(i-1)*3+2);
      patch([2 2 2.3 2.3], patchy, patchColor,'EdgeColor','none'); hold on;
      plot(tfit, octData, 'Color', name2rgb('red'), 'LineWidth', 1,'Clipping','off');
      xlim([tfit(1) tfit(end)]);
      ylim(yl);
      box on;
      set(gca,'xticklabel',[],'yticklabel',[]);
    end
    
    if (plotData(i,2))
      subplotp(Q,(i-1)*3+3);
      patch([2 2 2.3 2.3], patchy, patchColor,'EdgeColor','none'); hold on;
      plot(tfit, citData, 'Color', name2rgb('ForestGreen'), 'LineWidth', 1);    
      xlim([tfit(1) tfit(end)]);
      ylim(yl);
      box on;
      set(gca,'xticklabel',[],'yticklabel',[]);
    end
  end
  
%% Now annotate the plots
  titleFontSize  = 14;
  xlabelFontSize = 12;
  ylabelFontSize = 12;

  TITLE_  = @(x)  title(x,'FontSize',  titleFontSize);
  XLABEL_ = @(x) xlabel(x,'FontSize', xlabelFontSize);
  YLABEL_ = @(x) ylabel(x,'FontSize', ylabelFontSize);

  % Make the titles
  titlePlots = [1 4 7 10 13 16];
  titles     = {'Data', 'Constant Fit', 'Only Octanol', 'Only Citral*', 'Scaled Cit+Oct', 'Free'};  
  for i = 1:numel(titlePlots)
    subplotp(Q, titlePlots(i)); 
    TITLE_(titles{i});
  end
  
  % Label the x and y ticks on one of the axes
  subplotp(Q, 3); 
  xtick = (0:3)+2; xtickLabel = arrayfun(@num2str, xtick-2, 'UniformOutput', false);
  ytick = [0 1.5]; ytickLabel = arrayfun(@num2str, ytick,   'UniformOutput', false);
  set(gca,'xtick', xtick, 'xtickLabel', xtickLabel);
  set(gca,'ytick', ytick, 'ytickLabel', ytickLabel);
  XLABEL_('Time (s)'); 
  YLABEL_('FR/20 (Hz)');
  
  % Insert the model descriptions
  modelDescr = cellfun(@(x) strrep(x, 'y(t)', 'fit(t)'), modelDescr, 'UniformOutput', false);
  
  for i = 1:numel(modelDescr)
    modelDescr{i} = strrep(modelDescr{i},'x_1','oct');
    modelDescr{i} = strrep(modelDescr{i},'x_2','cit');
    indSplit = [1 strfind(modelDescr{i}, ' + ') numel(modelDescr{i})+1];
    descr    = {};
    if (i == 4)
      indSplit = [1 indSplit(2) indSplit(end)];
    end
    for j = 1:numel(indSplit)-1
      descr{j} = modelDescr{i}(indSplit(j):indSplit(j+1)-1);        
    end
    modelDescr{i} = descr;
  end  
  
  descrPlots = [4 8 12 15 18];
  for i = 1:numel(descrPlots)
    subplotp(Q, descrPlots(i)); 
    h = xlabel(modelDescr{i},'FontSize', 10);
    pos = get(h, 'Position');
    pos(1) = min(xlim);
    set(h,'HorizontalAlignment','left','Position',pos);
  end
end

%% Figure S2D: Model distribution
if (any(whichFigures==5))
  % Make pie charts showing the break down based on types of
  % response (Citral, Octanol, mixture).
  resultsFile = fullfile(figDir, thisDir, dataDir, 'resultsBestModel3Bayes1Laplace_slag100_sreg100_maxLag3_funcModelPrior.mat');
  Output = load(resultsFile);
  
  sfigure(FindFigureCreate('Figure S2D: Model distribution')); clf;
  set(gcf,'Color',[1 1 1],'Resize','off','NumberTitle','off');
  ResizeFigure(gcf,8,3,'inches');
  
  % The first element is the model types for the component, the second is the non-lagged models, the third is the lagged.
  modelTypes = {[2,9,5,12],... % octanal (unit, scaled), + lagged
                [3,10,6,13],... % citral  (unit, scaled), + lagged
                [4,11,7,14,8,15]}; % mixture  (unit, scaled, free), + lagged

  colors = cell(3,1); % component type, (unit/scaled/...).
  colors{1} = {'red3','red1','orangered','orange'}; % Unit, Unit-L, Scaled, Scaled-L
  colors{2} = {'green','green3','olivedrab3','olivedrab1'}; % Unit, Unit-L, Scaled, Scaled-L
  colors{3} = {'yellow3','yellow','gold1','lightgoldenrod2','skyblue3','skyblue1'};
  
  morphInds = 1:11;
  R = reshape(Output.Results, 168, []);
  R = R(:,morphInds);
  bm = cellfun(@(x) x.bestModel, R);
  bm = bm(:);
  responseTypeNames = {'oct','cit','mix'};
  responseTypeLabels = {'Octanol-Type','Citral-Type','Mixture-Type'};
  for i = 1:3    
    counts = [];
    for j = 1:numel(modelTypes{i})
      counts(j) = numel(find(bm==modelTypes{i}(j)));
    end
    if (i<3)
      fprintf('%s-type responses: %1.3e unit (total), %1.3e lagged (total)\n', responseTypeNames{i}, sum(counts(1:2))/sum(counts), sum(counts([2 4]))/sum(counts));
    end
    subplot(1,3,i);
    if (numel(counts)==4)
      labels = {'1','1-lag','k','k-lag'};
    else
      labels = {'1','1-lag','k','k-lag','f','f-lag'};
    end
    a = pie(counts, labels);
    cm = cellfun(@name2rgb, colors{i}, 'UniformOutput', false);
    cm = vertcat(cm{:});
    colormap(cm);
    freezeColors;
    text(0.5,-0.2,responseTypeLabels{i},'units','normalized','FontSize',14,'HorizontalAlignment','Center');
  end
  rgb2cm();
end

%% Figure S2E: Effect of mixture on fits
if (any(whichFigures==6))
  % Go through the results and grab, for the citral-type responses
  % the k values of the scaling, and plot them.
  resultsFile = fullfile(figDir, thisDir, dataDir, 'resultsBestModel3Bayes1Laplace_slag100_sreg100_maxLag3_funcModelPrior.mat');
  Output = load(resultsFile);

  morphInds = 1:11; % Indices in GetBinaryMixtureIndsForResponseLinearity for the binary mixture morphs

  R = reshape(Output.Results, 168, []);
  R = R(:,morphInds); % Grab just the data for the morphs.

  dataForFigure = load(fullfile(figDir, thisDir, dataDir, 'dataForFigure'));
  
  Icit     = dataForFigure.I{2}; % Index of cells best fit by citral response
  Ioct     = dataForFigure.I{3}; % Index of cells best fit by octanol response

  cellInds = {Ioct, Icit};
  
  unitModelInds   = {{2  9},{3 10}}; % Indices of the unlagged and lagged unit models, for octanol and citrial
  scaledModelInds = {{5 12},{6 13}}; % Indices of the unlagged and lagged scaled models.

  % Placeholders fo statistics used in the plots.
  mu    = zeros(2, numel(morphInds)-1); %mu(i,j) will hold the value for cellType i to mixture j.
  sd    = zeros(2, numel(morphInds)-1);
  se    = zeros(2, numel(morphInds)-1);
  pv    = zeros(2, numel(morphInds)-1);

  referenceMixture = [numel(morphInds) 1];
  mixtureOffset    = [0 1];

  kvals = {}; % Place holder to keep the scaling factors of the fits
  for cellType = 1:2
    I = cellInds{cellType};

    kvals{cellType} = nan(numel(I), numel(morphInds));

    % Go through each cell and get the scaling factor of its fits for
    % each mixture.
    
    for i = 1:numel(I)
      for j = 1:size(R,2)
        unitModels   = unitModelInds{cellType};
        scaledModels = scaledModelInds{cellType};
        switch(R{I(i),j}.bestModel)
         case unitModels   % Unit fit, unlagged or lagged
          kvals{cellType}(i,j) = 1; % Scaling factor is 1, by defn. of unit fit
         case scaledModels % Scaled Octanol, unlagged or lagged
          kvals{cellType}(i,j) = R{I(i),j}.w(2);
        end
      end
    end  
  
    % Loop through each of the mixtures and compute scaling effects
    for i = 1:numel(morphInds)-1
      % Find all the cells which had CELLTYPE response to both this
      % mixture and to the reference mixture. 
      ind      = find(~isnan(kvals{cellType}(:,referenceMixture(cellType))) & ~isnan(kvals{cellType}(:,i + mixtureOffset(cellType))));
      
      x1 = kvals{cellType}(ind, referenceMixture(cellType));  % Responses to oct140:cit30
      xn = kvals{cellType}(ind, i + mixtureOffset(cellType)); % Response  to this mixture.

      % Paired t-test to see if the means are different
      [h, pv(cellType, i)] = ttest(x1, xn); 
    
      mu(cellType, i)    = mean(xn-x1);
      sd(cellType, i)    = std(xn-x1);
      se(cellType, i)    = sd(cellType, i)/sqrt(numel(xn));
    end
  end
  
% Now plot everything 
  muOct = mu(1,:);
  muCit = mu(2,:);
  seOct = se(1,:);
  seCit = se(2,:);

  sfigure(FindFigureCreate('Figure S2E: Effect of mixture on fits')); clf;
  set(gcf,'Color',[1 1 1],'Resize','off', 'NumberTitle','off');
  ResizeFigure(gcf,8,3,'inches');

  % Plot citral
  plot(morphInds(1:end-1)-0.1,muCit,'o-','Color',name2rgb('ForestGreen'),'MarkerFaceColor',name2rgb('ForestGreen')); hold on;  
  line([1 1]'*[morphInds(1:end-1)]-0.1,[muCit-seCit;muCit+seCit],'Color',name2rgb('ForestGreen'));
  % Fit a sinusoid to the means
  u = sin(2*pi*(-1.2+(1:10))/9);
  w = [ones(10,1) u(:)]\muCit(:);
  muCite = [ones(10,1) u(:)]*w;
  plot(morphInds(1:end-1)-0.1, muCite,'--','Color',name2rgb('ForestGreen'));
 
  % Plot Octanol
  plot(morphInds(1:end-1)+0.1, muOct(end:-1:1),'o-','Color',name2rgb('red'),'MarkerFaceColor',name2rgb('red')); hold on;
  line([1 1]'*[morphInds(1:end-1)]+0.1,[muOct(end:-1:1)-seOct(end:-1:1);muOct(end:-1:1)+seOct(end:-1:1)],'Color',name2rgb('red'));
  u = sin(2*pi*(-1.5+(1:10))/9);
  muOct = muOct(end:-1:1);
  w = [ones(10,1) u(:)]\muOct(:);
  muOcte = [ones(10,1) u(:)]*w;
  plot(morphInds(1:end-1)+0.1, muOcte,'--','Color',name2rgb('red'));

  mixtureVals = GetBinaryMixturePairedConcentrations;
  mixtureInds = GetBinaryMixtureIndsForResponseLinearity;
  mixtureVals = mixtureVals(mixtureInds(morphInds),:);
  mixtureVals = mixtureVals(2:end, [2 1]);
 
  xticks = 1:10;
  xtickLabels = arrayfun(@(xt) when(mod(xt,2), sprintf('%d:%d', mixtureVals(xt,1), mixtureVals(xt,2)),''), xticks, 'UniformOutput', false);
  set(gca,'ytick',[-0.2:0.1:0.1],'xtick', xticks, 'xtickLabels', xtickLabels, 'TickLength',[0.01 0.01]/2, 'FontSize', 11);
  xlim([0 morphInds(end)]);
  ylim([-0.25 0.15]);
  ylabel('K(x) - K(140:30)', 'FontSize', 12);
  xlabel('Mixture (relative to response type)', 'FontSize', 12);
  HorizontalLine(gca, 0, 'k:');

  % Plot '*' where the p-value is < 0.05
  pvCit = pv(2,:);
  pvOct = fliplr(pv(1,:));
  
  arrayfun(@(i) plot(i,      0.02, '*','Color', (pvCit(i)<0.05)*name2rgb('ForestGreen') + (pvCit(i)>=0.05)*[1 1 1]), 1:numel(pvCit));
  arrayfun(@(i) plot(i+0.25, 0.02, '*','Color', (pvOct(i)<0.05)*name2rgb('red')         + (pvOct(i)>=0.05)*[1 1 1]), 1:numel(pvOct));

  % Plot the legend
  ylgnd = 0.12;
  line(0.0+[1.05 1.45]-0.5, [1 1]*ylgnd, 'LineWidth', 1, 'Color', name2rgb('ForestGreen'));
  text(0.0+1, ylgnd,'cit','FontSize',12);

  line(1.5+[1.05 1.45]-0.5, [1 1]*ylgnd, 'LineStyle', '--', 'LineWidth', 1, 'Color', name2rgb('ForestGreen'));
  text(1.5+1, ylgnd,'cit (fit)','FontSize',12,'HorizontalAlignment', 'left');

  line(3.0+[1.05 1.45]-0.5, [1 1]*ylgnd, 'LineStyle', '-', 'LineWidth', 1, 'Color', name2rgb('red'));
  text(3.0+1, ylgnd,'oct','FontSize',12,'HorizontalAlignment', 'left');

  line(4.5+[1.05 1.45]-0.5, [1 1]*ylgnd, 'LineStyle', '--', 'LineWidth', 1, 'Color', name2rgb('red'));
  text(4.5+1, ylgnd,'oct (fit)','FontSize',12,'HorizontalAlignment', 'left');
  
  TightenAxesToFigure;
end

%% Figure S2F: Response preference
if (any(whichFigures==7))
  % Plot the SNR angles for responses of each type.
  data = load(fullfile(figDir, thisDir, dataDir, 'snrAngleData'));
  th = data.th;
  ct = data.ct;
  
  colors = [name2rgb('red');name2rgb('forestgreen'); name2rgb('gold')];
  
  sfigure(FindFigureCreate('Figure S2F: Response preference')); clf; 
  set(gcf,'Color',[1 1 1],'Resize','off','NumberTitle','off');
  ResizeFigure(gcf,8, 4,'inches');
  
  Q = ComputeSubplotPosition(0.08,0.01,0.015,0.125);
  subplotp(Q,1);
  
  bw = [0.09 0.09 0.15];
  xi = -pi:pi/64:pi;
  h = [];
  for i = [3 1 2]    
    f = ksdensity(th(ct==i),xi,'width',bw(i));   
    h(i) = area(xi,f*first(diff(xi)));
    arrayfun(@(hh) set(hh,'FaceColor',colors(i,:),'EdgeColor',colors(i,:)), h(i));
    arrayfun(@(hh) set(get(hh,'Children'),'FaceAlpha',0.2), h(i));
    hold on;
  end
  legend(h([2 1 3]), 'cit-type','oct-type','mix-type');
  xlim([-pi/2 pi]);
  set(gca,'xtick',[-pi/2:pi/4:pi],'xticklabel',{'-pi/2','-pi/4','0','pi/4','pi/2','3pi/4','pi'},'ytick',[0:0.02:0.08],'FontSize',10);
  xlabel('SNR Angle','FontSize',12);
  ylabel('Frequency Observed','FontSize',12);
  ylim([-0.0004 0.08]);
  box on;
end