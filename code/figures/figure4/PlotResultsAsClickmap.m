function [colorMaps, colorMapNames, pnPlotOrder, I1, Imix, Iamb, Inull] = PlotResultsAsClickmap(resultsFile, pnCbotDataFileLong, pnCbotDataFile)
% [colorMaps, colorMapNames, pnPlotOrder, I1, Imix, Iamb, Inull] = PlotResultsAsClickmap(resultsFile, pnCbotDataFileLong, pnCbotDataFile)
%
% Plots the results of the model fitting procedure (contained in
% RESULTSFILE) as a 'click-map': CELLS x MIXTURES matrices colored
% according to the type of fit selected, or the R2 or SNR value of the
% fits. Clicking on a matrix element will produce a plot showing the
% corresponding fit, which is useful for manually checking the results
% of the procedure. PNCBOTDATAFILELONG is the name of the file
% containing the responses over a longer post odor offset period (used
% for computing SNR), and PNCBOTDATAFILE is the file containing the
% responses as used in the fits.
%
% The cells are grouped into four groups based on the non-constant
% fits to their responses, as follows:
%
% 1) Those for which no non-constant responses were found. Their
% indices are returned in INULL.
%
% 2) Those for which the majority of their non-constant responses
% required more than a single component to fit (IMIX). 
%
% 3) Those for which the majority of non-constant responses were fit
% by single component responses, but in which no single component
% dominated (IAMB).
%
% 4) Those for which the majority of the non-constant responses were
% fit by (the same) single component response. (I1).
%
% The index of the PN plotted in each row is returned in PNPLOTORDER.
%
% COLORMAPS is a cell array of matrices containing ((NUMCELLS x
% NUMMIXTURES) x 3) colormaps. The colors in the rows of each colormap
% index each element of the click-map, so that plotting the numbers
% 1:NUMCELLS x NUMMIXTURES as a NUMCELLS x NUMMIXTURES matrix and
% using a colormap will color each matrix according to the results for
% the corresponding cell and mixture.
%
% Four color maps are generated:
%
% 1) Colored by components used in the fit:
%    A: Red,  B: orange,  C: green,  D: brown,
%    W: blue, X: magenta, Y: yellow, Z: light green,
%    More than 1 component used: white, 
%    No components used (constant model): black.
%
% 3) Colored by R^2: Red: 0 -> 1: Blue
% 
% 2) Colored by component color (as in 1), scaled by R^2.
%
% 4) Colored by SNR: dark blue: low -> white: high.
%
% The names of the colormaps are returned in the cell array
% COLORMAPNAMES.
%
% Pressing space cycles through the colormaps, and clicking on a
% matrix element will show the corresponding fit.

ResultsStruct = load(resultsFile);
Results = ResultsStruct.Results;

fitOptions = {'sigmaLag', ResultsStruct.sigmaLag, ...
              'sigmaReg', ResultsStruct.sigmaReg, ...
              'lagLimit', ResultsStruct.lagLimit, ...
              'variancePriorAlpha', ResultsStruct.variancePriorAlpha,...
              'variancePriorBeta',  ResultsStruct.variancePriorBeta};

bestModel      = cellfun(@(x) x{1}, Results);
numCmpsInModel = cellfun(@(x) sum(x{2}), Results);

B  = GetOdorNamesAsBinaryVectors('full');
r2 = cellfun(@(x) x{3}(2), Results);

snr = ComputeSnrForAllCellsAndMixtures('dataFile', pnCbotDataFileLong);
snr = snr(:,1:size(Results,2));
snr(isinf(snr)) = max(snr(~isinf(snr)));
snrDb = 10*log10(snr); % Some might be -inf because snr = 0

colorMaps = {};
%% Make the color maps
% First sort the cells by their single component responses
cmps = 'ABCDWXYZ';
odors = GetOdorsList;
odors = cellfun(@(x) x(5:end),odors,'UniformOutput',false); odors{end} = cmps;
ic = cellfun(@(x) x{2}, Results,'UniformOutput',false);  % Grab the input configuration for each result
icn = zeros(size(Results));
% Go through each of the ICs and map them to a number indicating the
% components used in the fit.
for i = 1:size(Results,1) % cells
  for j = 1:size(Results,2) % mixtures
    switch(sum(ic{i,j}))
     case 0 % No components used: constant model.
      icn(i,j) = 0;
     case 1 % A single component used, record which one.
      thisCmps = odors{j+12};
      thisCmp = thisCmps(find(ic{i,j}));
      icn(i,j) = find(thisCmp==cmps);
     otherwise % More than one component used.
      icn(i,j) = -1;
    end
  end
end

[U,I] = OrderResponses(icn,[1:8 -1], 0);
pnPlotOrder = vertcat(I{:});
I1   = I(1:8);
Imix = I(9);
Iamb = I(10);
Inull= I(11);

pnGroupSizes = [numel(vertcat(I1{:})) numel([vertcat(Imix{:});vertcat(Iamb{:})]) numel(vertcat(Inull{:}))];
pnGroups = {1:pnGroupSizes(1), pnGroupSizes(1)+(1:pnGroupSizes(2)), pnGroupSizes(1)+pnGroupSizes(2)+(1:pnGroupSizes(3))};

% Colormap 1: Color by Single Component Response
componentColors = [name2rgb('white'); zeros(1,3); GetComponentColors(1)];
cols = zeros([size(U) 3]);
for i= 1:size(U,1)
  for j = 1:size(U,2)    
    cols(i,j,:) = componentColors(U(i,j)+2,:);
  end
end
colorMaps{1} = reshape(cols,[],3);
colorMapNames{1} = 'By Single Component Response';

% Colormap 2: Color by Single Component Response * R2
componentColors = [name2rgb('white'); zeros(1,3); GetComponentColors(1)];
cols = zeros([size(U) 3]);
if (isempty(r2))
  r2 = ones(size(U));
end
r2(r2<0) = 0;
for i= 1:size(U,1)
  for j = 1:size(U,2)    
    hsvCol = rgb2hsv(componentColors(U(i,j)+2,:));
    hsvCol(3) = hsvCol(3)*r2(pnPlotOrder(i),j);
    rgbCol    = hsv2rgb(hsvCol);
    cols(i,j,:) = rgbCol;
  end
end
colorMaps{2} = reshape(cols,[],3);
colorMapNames{2} = 'By Single Component Response x R2';

% Colormap 3: Color by Snr
cols = zeros([size(snrDb) 3]);
cols = interp1([-5; 0; 3; 5; 10; 15; 20; max(snrDb(:))],...
               [name2rgb('royalblue');...
                name2rgb('black');...                
                name2rgb('darkred');...
                name2rgb('red');...
                name2rgb('yellow');...
                name2rgb('white');...
                name2rgb('white');...
                name2rgb('white');],snrDb(:));
cols = reshape(cols,[size(U) 3]);
cols = cols(pnPlotOrder,:,:);
colorMaps{3} = reshape(cols,[],3);
colorMapNames{3} = 'By SNR';

% Colormap 4: Color by R2 
componentColors = [name2rgb('white'); zeros(1,3); GetComponentColors(1)];
cols = zeros([size(U) 3]);
r2(r2<0) = 0;

cols = interp1([0; 0.1; 0.25; 0.5; 0.75; 0.9; 1.0],...
               [name2rgb('black');...
                name2rgb('red');...
                name2rgb('orange');...
                name2rgb('white');...
                name2rgb('deepskyblue2');...
                name2rgb('royalblue3'); 
                name2rgb('royalblue3')],r2(:));
cols = reshape(cols,[size(U) 3]);
cols = cols(pnPlotOrder,:,:);
colorMaps{4} = reshape(cols,[],3);
colorMapNames{4} = 'By R2';

% Done with colormaps, now plot
maxColorModes = length(colorMaps);
colorMode = 4;
makeClickmap();
set(gcf,'KeypressFcn',@toggleColorMode);

function makeClickmap()
sfigure(FindFigureCreate('Complex Mixture Clickmap')); clf;
Mmax = prod(size(Results));
M = reshape(1:Mmax,size(Results));
Q = ComputeSubplotPositionsForMixturePlots(pnGroupSizes,32,0.01,0.8);
for i = 1:3
  if (pnGroupSizes(i)==0)
    continue;
  end
  subplotp(Q,i);
  set(imagesc([1 size(Results,2)], [pnGroups{i}(1) pnGroups{i}(end)], M(pnGroups{i},:),[1 Mmax]),'hittest','off');
  colormap(colorMaps{colorMode});
  set(gca,'ButtonDownFcn', @(h,evt) mapPointToCellMixture(get(gca,'CurrentPoint')));
  set(gca,'ytick',1:pnGroupSizes(i),'yticklabel',arrayfun(@num2str, pnPlotOrder(pnGroups{i}),'UniformOutput',false), 'FontSize',6,'ticklength',[0 0]);
  axis equal; axis tight;
end
set(gcf,'Name',sprintf('Colormap %d: %s', colorMode, colorMapNames{colorMode}));
end

function mapPointToCellMixture(pt)
whichCell = pnPlotOrder(round(pt(3)));
whichOdor = round(pt(1));
whichInputConfig = Results{whichCell, whichOdor}{2};
whichModel = Results{whichCell, whichOdor}{1};
ShowModelFits(whichCell, whichOdor+12, whichInputConfig, whichModel, pnCbotDataFile, 'snrDb', snrDb(whichCell, whichOdor), fitOptions{:});
end

function toggleColorMode(src, evt)
if (evt.Character == ' ')
  colorMode = colorMode+1;
  if (colorMode>maxColorModes)
    colorMode = 1;
  end
end
makeClickmap();
end

function str = ResultsString(lpost, sse, snrDb)
str = sprintf('LPOST = %1.2f   SSE = %1.3f   SNR = %1.3f dB', lpost, sse. snrDb);
end

end