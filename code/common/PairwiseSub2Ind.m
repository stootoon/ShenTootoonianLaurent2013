function ind = PairwiseSub2Ind(u,r,c)
% ind = PairwiseInd2Sub(u,r,c)
%
% Returns the index into a PDIST vector of the comparison between
% elements R and C, of a set of U items being compared.
if (r == c)
  error('The row and column indices cannot be equal.');
end
if (r>u | c>u)
  error('The row and column indices must be <= the number of items being compared.');
end
s = zeros(u);
s(r,c) = 1;
s(c,r) = 1;
ind = find(squareform(s));

