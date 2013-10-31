function [colorMap, I] =  PlotBayes1LaplaceResultsAsClickMap(Results, pnCbotDataFile, varargin)
%  [colorMap, I] = PlotResultsAsClickMap(Results, pnCbotDataFile, varargin)
%
% Given the NUMCELLS x NUMMIXTURES array of Results structures
% augmented with nonlinearity data (as returned by
% COMPUTELINEARITYINDICESFROMRESULTS), organizes and colorcodes the
% cell data and creates a click map showing the responses, and
% allowing the user to click on a cell and see its full response and
% fits pop up (using SHOWMODELFITSBAYES). Pressing any key will cause
% the plot to switch colorcodes from by-color to by-nonlinearity.
% 
% This function is based heavily on the script
% collectComputeBestModels2BayesResults.m

p = inputParser;
p.addOptional('pvals',[]);
p.addOptional('plotConcSeries',false);
p.addOptional('fitOptions',{});
p.addOptional('r2',[]);
p.addOptional('e2',[]);
p.addOptional('snrDb',[]);
p.parse(varargin{:});

r2 = p.Results.r2;
e2 = p.Results.e2;
snrDb = p.Results.snrDb;

fitOptions = p.Results.fitOptions;
pvals = p.Results.pvals;
plotConcSeries = p.Results.plotConcSeries;

[mixtureInds, mixtureVals] = GetBinaryMixtureIndsForResponseLinearity;
morphInds = 1:11;

numMixtures = numel(mixtureInds);
numCells = 168;

bestModel = cellfun(@(C) C.bestModel,Results);
bmHasOct  = cellfun(@(C) sum(C.indInputs==1)>0, Results);
bmHasCit  = cellfun(@(C) sum(C.indInputs==2)>0, Results);
octFlag   = (bmHasOct & ~bmHasCit);
citFlag   = (bmHasCit & ~bmHasOct);
mixFlag   = (bmHasCit & bmHasOct);
constFlag = (bestModel==1);

bestModel = reshape(bestModel,numCells, numMixtures);
citFlag   = reshape(citFlag,  numCells, numMixtures);
octFlag   = reshape(octFlag,  numCells, numMixtures);
mixFlag   = reshape(mixFlag,  numCells, numMixtures);
Results   = reshape(Results,  numCells, numMixtures);

nzSum  = sum(bestModel(:,morphInds)~=1,2);
zSum   = sum(bestModel(:,morphInds)==1,2);
citSum = sum(citFlag(:,morphInds),2);
octSum = sum(octFlag(:,morphInds),2);
mixSum = sum(mixFlag(:,morphInds),2);

s = [zSum citSum octSum mixSum];
cellType = zeros(size(s,1),1);
for i = 1:size(s,1)
  if (zSum(i) == numel(morphInds))
    cellType(i) = 1;
  else
    mx = argmax(s(i,2:end));
    if (length(mx)>1)
      cellType(i) = 5;
    else
      cellType(i) = mx+1;
    end
  end
end
%% Make the colormaps
colorModeValid = zeros(1,4);
colorMap = {};

%% Colormap 1: Based on best Model
colors = GetBayes1LaplaceModelColors(2);
bm = bestModel(:);
cmModels = zeros(size(bm,1),3);
for i = 1:15
  cmModels(bm==i,:) = ones(sum(bm==i),1)*colors(i,:);
end
colorMap{1} = cmModels;
colorModeValid(1) = 1;

%% Colormap 2: Based on coefficients
cmCoefs = 0*cmModels;
wmin = -1.5;
wmax = 1.5;
for i = 1:numCells
  for j = 1:numMixtures
    ind = (j-1)*numCells+i;
    col = [0 0 0];
    thisBm = bestModel(i,j);
    switch thisBm
     case 1 % Constant
      col = rand*0.1*[1 1 1];
     case {2,9}  % Unit Octanol
      col = [1 0 0];
     case {3,10}  % Unit Citral      
      col = [0 1 0];
     case {4,11}  % Unit Mixture
      col = [1 1 0];
     case {5,12}  % Scaled Octanol
      w = remapWeightIntoUnitInterval(Results{i,j}.w(2), wmin, wmax);
      col = [w 0 0];
     case {6,13}  % Scaled Citral
      w   = remapWeightIntoUnitInterval(Results{i,j}.w(2), wmin, wmax);
      col = [0 w 0];
     case {7,14}  % Scaled Mixture
      w = remapWeightIntoUnitInterval(Results{i,j}.w(2), wmin, wmax);
      col = w*[1 1 0];     
     case {8,15}  % Free Mixture
      w = arrayfun(@(u) remapWeightIntoUnitInterval(u, wmin, wmax), Results{i,j}.w(2:end));
      col = [w 0];     
    end
    cmCoefs(ind,:) = col;
    ind = ind + 1;
  end
end   
colorModeValid(2) = 1;
colorMap{2} = cmCoefs;

%% Colormap 3: Based on r2 and snrdB
if (~isempty(snrDb) && ~isempty(r2))
  redHue = 0; 
  blueHue= 2/3; 
  Mnl= cat(3, r2, snrDb);
  r2col = interp1([0 0.1 0.2 0.3 1],[1 0 0; name2rgb('red'); 1 1 1; 0 0 1; 0 0 1],r2);
  gain = arrayfun(@(x) mean(snrDb(:)<=x), snrDb(:));
  cmNl = bsxfun(@times, reshape(r2col,[],3), gain(:));
  colorModeValid(3) = 1;
  colorMap{3} = cmNl;
end

%% Colormap 4: Based on Best Single Component
cmComps = 0*cmNl;
bm = bestModel(:);
cmComps = (bm==1)*[0 0 0] +...
         (bm==2 | bm==5 | bm==9 | bm==12)*[1 0 0] +...
         (bm==3 | bm==10| bm==6 | bm==13)*name2rgb('ForestGreen')+...
         (bm==4 | bm==11| bm==7 | bm==14 | bm==8 | bm==15)*name2rgb('Gold');
colorModeValid(4) = 1;
colorMap{4} = cmComps;

%% Colormap 5: Based on R2
cmR2 = 0*cmNl;

cols = zeros([size(r2) 3]);
r2(r2<0) = 0;
cols = interp1([0; 0.1; 0.25; 0.5; 0.75; 0.9; 1.0],...
               [name2rgb('black');...
                name2rgb('red');...
                name2rgb('orange');...
                name2rgb('white');...
                name2rgb('deepskyblue2');...
                name2rgb('royalblue3'); 
                name2rgb('royalblue3')],r2(:));
colorMap{5} = reshape(cols,[], 3);
colorModeValid(5) = 1;

%% Colormap 6: Based on SNR
cmSnr = 0*cmNl;

cols = zeros([size(snrDb) 3]);
cols = interp1([min(snrDb(:));-10; -5; 0; 5; 10; 15; 20; 25; max(snrDb(:))],...
               [name2rgb('royalblue');...
                name2rgb('royalblue');...
                name2rgb('royalblue');...
                name2rgb('black');...
                name2rgb('red');...
                name2rgb('orange');...
                name2rgb('yellow');...
                name2rgb('white');...
                name2rgb('white');...
                name2rgb('white');],snrDb(:));
colorMap{6} = reshape(cols,[], 3);
colorModeValid(6) = 1;

%% Rearrange the colormaps to reflect the plot order
% First, compute the plotorder
I = cell(1,5);
counts = zeros(1,5);
for i = 1:5
  I{i} = [];
  inds =find(cellType == i);
  bm = bestModel(cellType == i,:);
  counts(i) = size(bm,1);
  if (isempty(bm))
    continue;
  end
  % Sort the cells by their center of response mass 
  centerOfMass = bsxfun(@rdivide, bm(:,morphInds)~=1, sum(bm(:,morphInds)~=1,2))*morphInds(:);
  [foo,p] = sort(centerOfMass);
  I{i} = inds(p);
end
pnPlotOrder = vertcat(I{:});
pnGroupSizes = cellfun(@numel, I);

% Now reorder the colormaps
for i = 1:numel(colorMap)
  if (colorModeValid(i))
    cm = reshape(colorMap{i},numCells,[],3);
    colorMap{i} = reshape(cm(pnPlotOrder,:,:),[],3);
  end
end

%% Now make the clickmap
plotMode.colorModeValid = colorModeValid;
plotMode.currentMode = 6;
plotMode.colorMaps = colorMap;
makeClickMap();

set(gcf,'ButtonDownFcn','', 'KeyPressFcn',@toggleColorModeAndPlot);

function toggleColorModeAndPlot(src, evt)
if (evt.Character == ' ')
  while(1)
    plotMode.currentMode = plotMode.currentMode + 1;
    if (plotMode.currentMode > length(plotMode.colorModeValid))
      plotMode.currentMode = 1;
    end
    if (plotMode.colorModeValid(plotMode.currentMode))
      break;
    end
  end
end
makeClickMap();
end

function makeClickMap()
fprintf('Plotting click map in "%d" mode.\n', plotMode.currentMode);

[numCells, numMixtures] = size(bestModel);

sfigure(FindFigureCreate('Binary Mixture Responses ClickMap')); clf;
set(gcf,'Resize','off');
ResizeFigure(gcf,11,8.5,'inches');
MAPPOINTTOCELLMIXTURE= @(pt,ind) ShowModelFitsBayes1Laplace(ind(round(pt(3))), round(pt(1)),...
                                                  bestModel(ind(round(pt(3))),round(pt(1))), pnCbotDataFile, fitOptions{:});

titles = {'Neither','Citral','Octanol','Mixture','Ambiguous'};
Q = ComputeSubplotPositions(1,5,[],0.05,0.05,0.05,0.05,0.075,0.0);
counts = 0*(1:5);
for i = 1:5
  subplotp(Q,i);
  M = reshape(1:prod(size(bestModel)),size(bestModel));
  hh = imagesc(M(sum(pnGroupSizes(1:i-1))+(1:pnGroupSizes(i)),1:when(plotConcSeries,16,11)),[1 numCells*numMixtures]);
  colormap(plotMode.colorMaps{plotMode.currentMode});
  freezeColors;

  set(gca,'ytick',1:pnGroupSizes(i),'yticklabel',arrayfun(@num2str,I{i},'UniformOutput',false),'FontSize',4);
  if (plotConcSeries)
    set(gca,'xtick',[1 6 11 12 14 16],'xticklabel',{'o30:c140','140:140','140:30','30:30','80:80','140:140'},'FontSize',6);
  else
    set(gca,'xtick',[1 6 11],'xticklabel',{'o30:c140','140:140','140:30'},'FontSize',6);    
  end
  xticklabel_rotate;
  set(gca,'ticklength',[0.005, 0.005]);
  set(hh,'HitTest','off');
  set(gca,'ButtonDownFcn',@(h,evt) MAPPOINTTOCELLMIXTURE(get(gca,'CurrentPoint'), I{i}));
  set(gca,'FontSize',6);
  set(gca,'Units','normalized');
  axis equal;
  axis tight;
  title(sprintf('%s (%d PNs)',titles{i},pnGroupSizes(i)),'FontSize',8);
end
end

function w1 = remapWeightIntoUnitInterval(w, minw, maxw)
w1 = (w-minw)/(maxw-minw);
w1 = min(w1,1);
w1 = max(w1,0);
end
end