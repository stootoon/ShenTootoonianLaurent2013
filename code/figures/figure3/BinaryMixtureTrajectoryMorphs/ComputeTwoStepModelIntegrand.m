function q = ComputeTwoStepModelIntegrand(th1, th2, x, y, dy)
% q = ComputeTwoStepModelIntegrand(th1, th2, x, y, dy)
%
% For given threshold levels th1 and th2, computes the integral over
% y1, y2, y3, and v of the y and v-dependant part of the integrand for the
% posterior of the two step model.
%
% The y and v dependent integrand is 
%
% v^(-k/2-1);
%
% where the dependence on y comes in through v. This is the likelihood
% of the parameters, i.e. P(D|theta,M), without the constant factors.

if (th1>=th2)
  q = 0;
  return;
end

i1 = x < th1;
i2 = (x < th2) & ~i1;
i3 = ~(i1 | i2);

k = numel(x);

k1 = sum(i1);
k2 = sum(i2);
k3 = sum(i3);

% Must be a non-degenerate two-step function
if (~k1 | ~k2 | ~k3) 
  q = 0;
  return;
end

y1 = mean(y(i1));
y2 = mean(y(i2));
y3 = mean(y(i3));

r1 = sum((y(i1)-y1).^2);
r2 = sum((y(i2)-y2).^2);
r3 = sum((y(i3)-y3).^2);

v = (r1+r2+r3)/(k+2);

f_pk = v^(-k/2-1)*exp(-(r1+r2+r3)/(2*v));

Hdiag = -1/v*[k1 k2 k3 (2+k)/(2*v)]; % diagonal of the hessian. Off-diagonal elements are zero.
iHdiagNz = Hdiag~=0;
nd = sum(iHdiagNz);
volGauss = (sqrt(2*pi)^nd)/sqrt(prod(abs(Hdiag(iHdiagNz))));
volCube = dy^(4-nd); % The diagonal term for the variance will never be zero, so the length term will always be dy.
q = f_pk * volGauss * volCube;

