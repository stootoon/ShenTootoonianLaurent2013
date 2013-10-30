function cols = GetMixtureResponseLinearityColors(whichMode)
% cols = GetMixtureResponseLinearityColors(whichMode)
%
% Returns the colors used for making the colormaps etc. for the
% mixture response linearity figures.

switch(whichMode)
 case 1 % By best model
  cols = [name2rgb('gray10');
          name2rgb('blue');
          name2rgb('magenta');
          name2rgb('red');
          name2rgb('deepskyblue3');
          name2rgb('violet');
          name2rgb('salmon')];
 case 2 % By num components
  cols = [name2rgb('gray10');
          name2rgb('cobalt'); % 1
          name2rgb('emeraldgreen'); % 2
          name2rgb('gold');     % 3
          name2rgb('orange'); % 4
          name2rgb('red'); % 5
          name2rgb('darkred'); % 6
          name2rgb('orchid4'); % 7 
          name2rgb('orchid1')]; % 8
 case 3 % By component
  cols = GetComponentColors;
 otherwise
  error('Unknown mode: %d', whichMode);
end
