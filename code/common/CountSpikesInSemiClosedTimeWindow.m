function F = CountSpikesInSemiClosedTimeWindow(spt,t0,t1)
% F = CountSpikesInSemiClosedTimeWindow(spt,t0,t1)
% 
% Accepts a matrix containing spike times in its columns and returns it 
% as an integer matrix containing the number of spikes in the [t0 t1) window.
%
% If SPT is sparse, it's made full and all entries set to zero
% are set to -Inf.
%
% Example:
%
% spt = [1 2 3; -1 2 3]'
%
% spt =
%      1    -1
%      2     2
%      3     3
% >> CountSpikesInSemiClosedTimeWindow(spt,0,2)
%
% ans =
%
%      1     0


IP = inputParser;
IP.addRequired('spt',   @(x) validateattributes(x, {'numeric'}, {'2d','nonnan','nonempty'}));
IP.addRequired('t0',    @(x) validateattributes(x, {'numeric'}, {'scalar','nonnan','nonempty'}));
IP.addRequired('t1',    @(x) validateattributes(x, {'numeric'}, {'scalar','nonnan','nonempty','>',t0}));
IP.parse(spt, t0, t1);

if (issparse(spt))
  spt = full(spt);
  spt(spt==0) = -Inf;
end

F = zeros(size(spt));
F(spt>=t0 & spt<t1) = 1;
F = sum(F,1);

