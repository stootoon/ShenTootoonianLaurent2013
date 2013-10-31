function str = GetCodeDirForFigure(whichFigure)
% str = GetCodeDirForFigure(whichFigure)
%
% Returns a string containing the code directory for the specified
% figure number. WHICHFIGURE must be in the range 2-8.

if (whichFigure<2 | whichFigure>8)
  error('Figure must be in the range 2-8.');
end

str = fullfile(GetRootDir('figures'),['figure' num2str(whichFigure)]);
