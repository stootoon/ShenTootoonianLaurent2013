function cols = GetComponentColors(whichMode)
% cols = GetComponentColors(whichMode)
%
% Returns an 8 x 3 matrix containing RGB colors corresponding to the
% odors used in the mixture experiments.

if (nargin==0)
  whichMode = 1;
end
switch(whichMode)
 case 1 % Kai's style
  cols = [name2rgb('red');
          name2rgb('orange');
          name2rgb('ForestGreen');
          name2rgb('sienna');
          name2rgb('royalblue');
          name2rgb('mediumorchid2');
          name2rgb('gold');
          name2rgb('limegreen');];
 case 2
  cols = [name2rgb('darkred'); % Octanol: A
          name2rgb('red'); % 2-Nonanone: Y
          name2rgb('orange');  % Citral: C
          name2rgb('gold'); % Isoamyl cetat: W
          name2rgb('LimeGreen'); % Butanedione X
          name2rgb('ForestGreen'); % Phentole B
          name2rgb('cobalt');% Benzaldehyde D
          name2rgb('violet')]; % L-carvone Z
  [foo, ord] = sort([1 7 3 5 6 2 4 8]);
  cols = cols(ord,:);
 case 3
  h = linspace(0,1,9);
  h = h([1:3 5:8 4]);
  colsHsv = [h',ones(8,1) ones(8,1)];
  
  cols = hsv2rgb(colsHsv);
end
