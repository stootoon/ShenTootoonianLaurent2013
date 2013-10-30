function [r,c] = PairwiseInd2Sub(u,ind)
% [r,c] = PairwiseInd2Sub(u,ind)
%
% Returns the row/column subscript values of the index IND into the
% vector of pairwise comparisons U. U for example is the result of a
% call to PDIST, and R and C are the row and column indices of the
% corresponding element after a call to SQUAREFORM. 
%
% If U is a scalar, it's treated as the the number of elements that
% were pairwise compared.

if (isscalar(u))
  u = pdist(ones(u));
else
  u = 0*u;
end

u(ind) = 1;
sq = squareform(u);
ind = find(sq);
[c,r] = ind2sub(size(sq), ind);
ind = find(c>r);
c = c(ind);
r = r(ind);

