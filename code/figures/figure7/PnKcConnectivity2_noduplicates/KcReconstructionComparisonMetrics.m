function m = KcReconstructionComparisonMetrics(X,Xrec,whichMetric,summaryFunction,varargin)
% m = KcReconstructionComparisonMetrics(X,Xrec,whichMetric,summaryFunction,varargin)
%
% This function is called by MAKEFIGURES and computes metrics for
% comparing KC responses to their reconstructions.
%
% X and Xrec should be NUMSAMPLES x NUMKCS matrices of KC responses
% and their reconstructions.
%
% whichMetric can be one of
%
% 'SSE': The sum of the squared error between the KC response
% and its reconstruction,
%
% 'R2': 1 minus the ratio of SSE/SST, where SST is the sum of squares
% of the KC response, 
%
% 'CorrDist': The correlation distance between the KC response
% and its reconstruction 
%
% 'EucDist': The euclidean distance between the KC response
% and its reconstruction,
%
% 'pvalue': The pvalue for the regression. An optional argument 'df'
% should provide a 1 x NUMKCS vector of degrees of freedom for the
% regression.
%
% Each of these functions produces a 1 x NUMKCS vector of metrics,
% computed for each of KCs. SUMMARYFUNCTION is applied to this matrix
% to produce a final summary metric. A typical value for it is
%
% @(m) mean(m).

p = inputParser;
p.addOptional('df',[]);
p.parse(varargin{:});

df = p.Results.df;

if (length(size(X))~=2)
  error('X should have two dimensions.');
end

if (size(X)~=size(Xrec))
  error('X and Xrec should have the same size.');
end



[numSamples,numCells] = size(X);

switch (lower(whichMetric))
 case 'sse'
  M = sum((X - Xrec).^2,1);
 case 'ssen'
  M = sum((X - Xrec).^2,1)./sum(bsxfun(@minus,X, mean(X)).^2);
 case 'r2'
  sst = sum(bsxfun(@minus, X, mean(X,1)).^2);
  sse = sum((X - Xrec).^2,1);
  M = 1 - sse./sst;
 case 'corr'
  Xm = bsxfun(@minus, X, mean(X,1));
  Xrecm = bsxfun(@minus,Xrec, mean(X,1));
  M = sum(Xm.*Xrecm)./normc(Xm)./normc(Xrecm);
 case 'modifiedr2'
  Xm = bsxfun(@minus, X, mean(X,1));
  Xrecm = bsxfun(@minus,Xrec, mean(X,1));
  Res = Xm - Xrecm;
  Nxm = normc(Xm);
  Nxrecm = normc(Xrecm);
  Nres = normc(Res);
  cosTheta = (Nres.^2 - Nxm.^2 - Nxrecm.^2)./(2*Nxm.*Nxrecm);
  M = cosTheta.^2;
 case 'corrdist'
  M = corrdist(X,Xrec);
 case 'eucdist'  
  M = normc(X-Xrec);
 case 'pvalue'
  if (isempty(df))
    error('No degrees of freedom supplied for p-value computation.');
  elseif numel(df)~= numCells
    error('df vector must have the same number of elements as cells.');
  elseif any(df==0)
    error('Degrees of freedom must be positive.');
  end
  mx = mean(X,1);
  ssr = sum(bsxfun(@minus, Xrec, mx).^2,1);
  sse = sum(bsxfun(@minus, Xrec, X).^2,1); 
  dfr = df - 1; 
  dfe = numSamples - df; 
  f = (ssr./dfr)./(sse./dfe); 
  M(dfr>0) = 1-fcdf(f(dfr>0),dfr(dfr>0),dfe(dfr>0));
  M(dfr==0) = 1;
 otherwise
  error('Unknown metric "%s".\n', whichMetric);
end

m = summaryFunction(M);

