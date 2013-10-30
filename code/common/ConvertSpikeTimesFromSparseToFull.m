function Y = ConvertSpikeTimesFromSparseToFull(X)
% Y = ConvertSpikeTimesFromSparseToFull(X)
%
% Converts the matrix of spiketimes X from sparse to full
% by setting all values that are 0 to -Inf.

IP = inputParser;
IP.addRequired('X', @(x) validateattributes(x, {'numeric'}, {'2d','nonempty','nonnan'}));
IP.parse(X);

Y = full(X); 
Y(Y == 0) = Inf;
Y(isinf(Y)) = Inf;
Y = sort(Y);
Y(isinf(Y)) = -Inf;

