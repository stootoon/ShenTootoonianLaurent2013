function [m, ci] = ComputeMedianAndConfidenceIntervalsFromData(X, alpha, pc)
% [m, ci] = ComputeMedianAndConfidenceIntervalsFromData(X, alpha, pc)
%
% Computes the medians for each of the columns of X, as well as the
% (1-alpha) confidence intervals, the latter using the large sample
% approxiation. If pc,the percentile is not provided, it defaults to
% 50. ALPHA should be in [0 1], PC should be [0 to 100].

if (nargin==2)
  pc = 0.5;
end

Y = X;
if (IsVector(X))
  Y = X(:);
else
  Y = reshape(X,size(X,1),[]);
end

m  = zeros(size(Y,2),1);
ci = zeros(size(Y,2),2);
z  = sqrt(2)*erfinv(alpha);
p  = pc/100;
for i = 1:size(Y,2)
  y = Y(:,i);
  y = y(~isnan(y));
  n = numel(y);  
  
  if (n)
    [foo,I] = sort(y);
    [foo,r] = sort(I);
    
    rm = n*p
    ru = (n*p + z*sqrt(n*p*(1-p)));
    rl = (n*p - z*sqrt(n*p*(1-p)));
    m(i)    = y(find(r==round(rm)));
    ci(i,1) = y(find(r==round(rl)));
    ci(i,2) = y(find(r==round(ru)));
  else
    m(i) = 0;
    ci(i,:) = [0 0];
  end
end
  

