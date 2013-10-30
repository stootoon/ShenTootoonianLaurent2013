function cols = GetBayes1LaplaceModelColors(whichMode)
% cols = GetBayes1LaplaceModelColors(whichMode)
%
% Returns colors for the various fit models. 

if (nargin == 0)
  whichMode = 1;
end

switch(whichMode)
case 1
 unlaggedColors = [name2rgb('red');     % Unit Octanol
                   name2rgb('green2');  % Unit Citral
                   name2rgb('yellow');   % Unit Mixture
                   name2rgb('crimson'); % Scaled Octanol
                   name2rgb('green3');  % Scaled Citral
                   name2rgb('gold');   % Scaled Mixture
                   1 1 0.6];  % Free Mixture

 cols = [0 0 0;          % Constant model
         unlaggedColors*0.9;
         unlaggedColors*0.6];
case 2
 cols           = [name2rgb('black');       % 1. Constant
                   name2rgb('darkred');         % 2. Unit Octanol
                   name2rgb('darkgreen');      % 3. Unit Citral
                   name2rgb('royalblue4');        % 4. Unit Mixture
                   name2rgb('red');         % 5. Scaled Octanol
                   name2rgb('green4');      % 6. Scaled Citral
                   name2rgb('royalblue3');       % 7. Scaled Mixture
                   name2rgb('magenta');      % 8. Free Mixture
                   name2rgb('orangered');      % 9. Unit Octanol
                   name2rgb('limegreen');  % 10. Unit Citral
                   name2rgb('deepskyblue3');       % 11. Unit Mixture
                   name2rgb('orange');        % 12. Scaled Octanol
                   name2rgb('lawngreen');  % 13. Scaled Citral
                   name2rgb('skyblue');       % 14. Scaled Mixture
                   name2rgb('violet')];    % 15. Free Mixture

 otherwise
  error('Unknown color mode %d.', whichMode);
end
