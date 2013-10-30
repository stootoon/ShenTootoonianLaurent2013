function t1 = CheckStartTimeAndSetEndTimeUsingResponseWindow(t0, responseWindow, binSize)
% t1 = CheckStartTimeAndSetEndTimeUsingResponseWindow(t0, responseWindow, binSize)
%
% This function checks to makes sure that t0 is <= responseWindow(1),
% its start, and sets t1 to responseWindow(2)+2*binSize. These
% parameters are used to set the start and end times of the Parameters
% structure for trajectory computations.

if (t0>responseWindow(1))
  error('T0 = %1.3f is greater than the start of the response window, %1.3f.', t0, responseWindow(1));
end
t1 = responseWindow(2)+2*binSize;
