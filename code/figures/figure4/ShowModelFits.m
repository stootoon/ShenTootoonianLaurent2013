function ShowModelFits(whichCell, whichOdor, whichInputConfig, whichModel, pnCbotDataFile, varargin)
% ShowModelFits(whichCell, whichOdor, whichInputConfig, whichModel, varargin)
%
% Fits the response of WHICHCELL to WHICHODOR using the model indexed
% by WHICHMODEL and the input configuration
% WHICHINPUTCONFIGURATION. Which odor should be an integer in the
% standard indexing (1 = odorAhigh), e.g. as returned by
% ODOR_NAME_TO_INDEX. WHICHINPUTCONFIGURATION and WHICHMODEL are a
% binary vector and an integer, respectively, as described in
% FITALLMODELSFOROBSERVATION2LAPLACEGENERAL. PNCBOTDATAFILE is the
% name of the file containing the odor responses. Optional arguments
% are passed on to FITALLMODELSFOROBSERVATION2LAPLACEGENERAL.
%
% See also: FITALLMODELSFOROBSERVATION2LAPLACEGENERAL

fitOptions = inputParser;
fitOptions.KeepUnmatched = true;
fitOptions.parse(varargin{:});

fitOptions = UnpackStructureFieldsAsNameValuePairs(fitOptions.Results);

Data = load(pnCbotDataFile);
whichBinsToFit = find(Data.iall==Data.ifit(1))+(0:numel(Data.ifit)-1);

%% Fit the response
[odorCmps, odorNames] = GetOdorNamesAsBinaryVectors('full'); 
thisOdorCmps = find(odorCmps(whichOdor,:)); % Find the components in this mixture
cmps = 'ABCDWXYZ';

Xobs = squeeze(Data.pnCbot(whichCell, :, thisOdorCmps+3, :)); % Grab the component responses

if (numel(thisOdorCmps)==1)
  % In case we're fitting a single component 'mixture', the squeeze
  % operation above will have removed the odor dimension, so
  % reintroduce it.
  Xobs = reshape(Xobs, size(Xobs,1), 1, size(Xobs,2)); 
end
yobs = squeeze(Data.pnCbot(whichCell, :, whichOdor, :)); % Grab the mixture response

% Now perform the fit.
Fit = FitAllModelsForObservations2LaplaceGeneral(Xobs, yobs, 'whichBinsToFit', whichBinsToFit, fitOptions{:}, 'fitMode', 'single', 'whichInputConfig',whichInputConfig,'whichModel',whichModel);

%% Plot the best fit
sfigure(FindFigureCreate('Complex Mixture Fits: Best Fit'));
figName = sprintf('Complex Mixture Fits: Best Fit for Cell %d, Mixture %d (%s)', whichCell, whichOdor, odorNames{whichOdor});
set(gcf,'name', figName);
clf;
Q = ComputeSubplotPositions(size(Xobs,2)+1,1,[],0.05,0.05,0.05,0.1,0.05,0.01);
subplotp(Q,1);

U     = [mean(Xobs,3) mean(yobs,2)]; % Trial averages
ylims = [-0.1 2.5];

plot(Data.tall, mean(yobs,2), 'k','LineWidth',2);
hold on;
ye  = Fit.Model.test(Fit.Model.X);
col = GetMixtureResponseLinearityColors(1);
h   = plot(Data.tfit, ye,'Color', name2rgb('blue'), 'LineWidth',2);
set(gca,'xticklabel',[]);
ylim(ylims);

descr = Fit.descr;
for i = 1:numel(Fit.inputConfigs)
  descr = strrep(descr,['x_' num2str(i)], cmps(thisOdorCmps(i))); 
end
titleStr = sprintf('Fit: %s\nR^2 = %1.3f, SSE = %1.3f, Lpost = %1.3f', descr, Fit.Model.r2, Fit.Model.sse, Fit.logPosteriors);

title(titleStr,'FontSize',12);
cmpColors = GetComponentColors;
for i = 1:numel(thisOdorCmps)
  col = cmpColors(thisOdorCmps(i),:);
  subplotp(Q,i+1);
  plot(Data.tall, U(:,i),'Color',col,'LineWidth',2);
  ylim(ylims);
  if (whichInputConfig(i))
    set(gca,'Color',name2rgb('lemonchiffon'));
  end
  if (i~=numel(thisOdorCmps))
    set(gca,'xticklabel',[]);
  end
end

drawnow;
