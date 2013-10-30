function [Y,I] = OrderResponses(X, respVals, noRespVal)
% [Y,I] = OrderResponses(X, respVals, noRespVal)
%
% Given an integer matrix of responses X, in which Xij is the response
% of the i'th cell to the j'th mixture, reorders the rows of X to
% group cells by the component response (respVal) that dominated the
% fits. For those models in which a single component response
% domimanted, it also sorts them by their center of mass with respect
% to this response. Center of mass is defined as the average of the
% column (mixture) indices in which the dominant fit type was used.
%
% RESPVALS is a vector containing the integers indicating responses of
% different types, and NORESPVAL is the scalar value that indicates no
% response.
%
% I is a cell array whose elements contain the indices of the cells in
% which each of the response types dominated, and then the cells in
% which no response type dominated, and finally the cells in which no
% responses were generated.
%
% Y is X but with rows arranged according to I.

[numCells, numMixtures] = size(X);
Xs = zeros(numCells, numel(respVals)+1);
for i = 1:numel(respVals)
  Xs(:,i) = sum(X==respVals(i),2);
end
Xs(:,end) = sum(X==noRespVal,2);
flagNoResp = Xs(:,end)==numMixtures; % # noRespVal == numMixtures => all responses were no response. 

mixtureRespVal = numel(respVals)+1;
noRespVal      = numel(respVals)+2; % Note that we're overwriting this value.

% Now go through and determine which response value dominanted each
% response.
maxVals = zeros(numCells,1);
for i = 1:numCells
  if (flagNoResp(i))
    maxVals(i) = noRespVal;
  else
    imx = argmax(Xs(i,1:end-1));
    if (numel(imx)==1) % A single response value dominated, so record which one.
      maxVals(i) = imx;
    else
      maxVals(i) = mixtureRespVal; % More than one response value dominated.
    end
  end
end

% Group the responses by the response type that dominated.
I = arrayfun(@(i) find(maxVals==i),1:numel(respVals)+2,'UniformOutput',false);

% Now loop through the single component responses and sort by center of mass
for i = 1:numel(respVals)
  ir = I{i}; % Indices of the cells using this response type
  U = double(X(ir,:)==respVals(i)); % Flag locations where the 
  Us = sum(U,2);
  pos = 1:size(U,2);
  centerOfMass = bsxfun(@rdivide, U, Us)*pos(:); % A weighted sum of the columns using that fit type.
  [foo,ii] = sort(centerOfMass);
  I{i} = I{i}(ii); % Rearrange the cells by their center of mass.
end

% Now return X but rearranged according to I.
iy = 1;
for i = 1:numel(respVals)+2
  ir = I{i};
  Y(iy+(0:numel(ir)-1),:) = X(ir,:);
  iy = iy+numel(ir);
end
