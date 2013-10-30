function [logp, D] = ComputePosteriorOfTwoStepModel(x,y,xr,yr,nmc)
% [logp, D] = ComputePosteriorOfTwoStepModel(x,y,xr,yr,nmc)
%
% Similar to COMPUTEPOSTERIOROFSINGLESTEPMODELDIRECTLY, but instead
% assumes a model where the first step is followed by a second step.
% The logarithm of the posterior probability is returned in LOGP, and
% information about the maximum likelihood fit, useful for plotting,
% is returned in the stucture D.

k = numel(x);
dx = x(end)-x(1); % Use the actual x range rather than xr to get the range of X values.
dy = diff(yr);

C = (2*pi)^(-k/2)/dx/dy^3; % Constant term for the integration

% likelihood function
flklhd = @(Th1, Th2) arrayfun(@(th1, th2) ComputeTwoStepModelIntegrand(th1, th2,x,y,dy), Th1, Th2);

% Sample uniformly over the integration region. 
Th = rand(nmc,2)*dx+xr(1); % we don't filter for th1<th2, the likelihood function takes care of this.
lklhd = flklhd(Th(:,1),Th(:,2));
q     = lklhd./(xr(end)-Th(:,1));
logp = log(C*mean(q(:))/dx/dx); % the additional dx^2 term is from the MC sampling.

% Find the approximate ML parameters for the model, which is where
% sigma2 is minimized. We can't do this using fminsearch because the
% error surface is block due to the spacing between datapoints, so
% we'll just use the best results from the montecarlo runs.

if (nargout == 2) 
  ii = find(lklhd==max(lklhd(:)),1);
  
  th1_pk = Th(ii,1);
  th2_pk = Th(ii,2);
  
  y1_pk = mean(y(x<th1_pk));
  y2_pk = mean(y(x>=th1_pk & x<th2_pk));
  y3_pk = mean(y(x >= th2_pk));

  D.th1_pk = th1_pk;
  D.th2_pk = th2_pk;
  
  D.y1_pk = y1_pk;
  D.y2_pk = y2_pk;
  D.y3_pk = y3_pk;
  
  D.fval = @(u) arrayfun(@(x) when(x<th1_pk,y1_pk,x>=th2_pk,y3_pk, y2_pk), u);
  
  D.v_pk = sum((D.fval(x)-y).^2)/(k+2);
end

