function ShowModelFitsBayes1Laplace(whichCell, whichMixture, whichModel, pnCbotDataFile, varargin)
% ShowModelFitsBayes1Laplace(whichCell, whichMixture, whichModel, pnCbotDataFile, varargin)

p = inputParser;
p.KeepUnmatched = true;
p.parse(varargin{:});

argsForFit = UnpackStructureFieldsAsNameValuePairs(p.Unmatched);

Data = load(pnCbotDataFile);
tall = Data.binStarts(Data.whichBinsAll);
tfit = Data.binStarts(Data.whichBinsToFit);

[mixtureInds, mixtureVals] = GetBinaryMixtureIndsForResponseLinearity();

Fit = FitAllModelsForCellAndMixtureBayes1Laplace(whichCell, whichMixture, 'pnCbotDataFile', pnCbotDataFile, argsForFit{:});

bestModel = Fit.bestModel;
if (bestModel ~= whichModel)
  disp('WARNING: The best model specified in input and best model in refit don''t match. Were the same fit parameters used?');
end

disp('Constant Model vs Best Model: ');
disp(ResultsString(Fit, 1));
disp(ResultsString(Fit, bestModel));

colors = GetBayes1LaplaceModelColors();

sfigure(FindFigureCreate('Model Fits')); clf; set(gcf,'Color',[0.8 0.8 0.8],'name',sprintf('Cell %d, Mixture %d:%d (%d)', whichCell, mixtureVals(whichMixture,1), mixtureVals(whichMixture,2), whichMixture));

Q = ComputeSubplotPositions(14,2,{{1:2,1:2},{3:4,1:2},{5:6,1:2},{7,1},{8,1},{9,1},{10,1},{11,1},{12,1},{13,1},{14,1},{7,2},{8,2},{9,2},{10,2},{11,2},{12,2},{13,2},{14,2}},0.025,0.01,0.01,0.01,0.01,0.01);
subplotp(Q,1);
numAllBins = size(Fit.Xall,1);

ymax = max(Fit.X(:));
ymax = max(ymax,0.5);

rawTraceColor = name2rgb('gray75');
imagesc([tall(1) tall(end)], [0 ymax],squeeze(Fit.Uall(:,:,1)'), [0 5]); hold on;
colormap(1-gray(100)); axis xy;
set(line([2 2; 2.3 2.3]',[-0.1 ymax]'*[1 1]),'Color',name2rgb('black'),'LineWidth',2);
plot(tall,Fit.Xall(:,1),'r','LineWidth',2); ylim([-0.1 ymax]); grid on;
set(gca,'xticklabel',[]);

subplotp(Q,2);
imagesc([tall(1) tall(end)], [0 ymax],squeeze(Fit.Uall(:,:,2)'),[0 5]); hold on;
colormap(1-gray(100)); axis xy;
set(line([2 2; 2.3 2.3]',[-0.1 ymax]'*[1 1]),'Color',name2rgb('black'),'LineWidth',2);
plot(tall,Fit.Xall(:,2),'Color',name2rgb('ForestGreen'), 'LineWidth',2); ylim([-0.1 ymax]); grid on;
set(gca,'xticklabel',[]);

subplotp(Q,3);
imagesc([tall(1) tall(end)], [0 ymax],squeeze(Fit.Uall(:,:,3)'),[0 5]); hold on;
colormap(1-gray(100)); axis xy;
set(line([2 2; 2.3 2.3]',[-0.1 ymax]'*[1 1]),'Color',name2rgb('black'),'LineWidth',2);
hdata = plot(tfit, Fit.X(:,3),'Color','b','LineWidth',2); hold on;
hfit = plot(tfit, Fit.Models{bestModel}.test(Fit.Args{bestModel}),'Color',colors(bestModel,:),'LineWidth',2); ylim([-0.1 ymax]);
xlim([tall(1) tall(end)]);
h = text(0.05,1,ResultsString(Fit, bestModel),'FontSize',14,'FontWeight','bold');
ph = get(h,'Position');
ph(2) = max(get(gca,'ylim'));
set(h,'Position',ph,'VerticalAlignment','top');
grid on;
set(legend([hdata, hfit], 'data',Fit.Models{bestModel}.descr),'box','off','color','none');
set(gca,'xticklabel',[]);

for i = 1:numel(Fit.Models) 
  subplotp(Q,3+(when(i<9,i,i+1)));
  patch([2 2.3 2.3 2],[-0.1 -0.1 ymax ymax],name2rgb('lavender'),'EdgeColor','none'); hold on;
  plot(tfit, Fit.X(:,3),'Color','b','LineWidth',2); 
  h = plot(tfit, Fit.Models{i}.test(Fit.Args{i}),'Color',colors(i,:),'LineWidth',2);
  set(legend(h, Fit.Models{i}.descr), 'box', 'off', 'color', 'none');
  ylim([-0.1 ymax]);
  grid off; box on;
  line([2 2.3; 2 2.3]', [-0.1 -0.1; ymax ymax]');
  
  h = text(first(xlim)+diff(xlim)*0.01,1,ResultsString(Fit, i),'FontSize',8);
  ph = get(h,'Position');
  ph(2) = max(get(gca,'ylim'));
  set(h,'Position',ph,'VerticalAlignment','top');
  %  title(ResultsString(Fit,i),'FontSize',8);
  set(gca,'xticklabel',[],'yticklabel',[],'ticklength',[0 0]);
  if (i==bestModel)
    set(gca,'Color',name2rgb('lemonChiffon'));
  end
end

function lstr = ResultsString(Fit, whichModel)
lstr = sprintf('Lpost = %1.3f SSE = %1.3f R^2 = %1.3f rho = %1.3f',Fit.logModelPosteriors(whichModel), Fit.sse(whichModel), Fit.r2(whichModel), Fit.rho(whichModel));
if (whichModel == Fit.bestModel)
  %lstr = sprintf('%s R^2 = %1.3f, SNR = %1.3f dB',lstr, result.r2, result.snrDb);
end