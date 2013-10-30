function str = GetDataDirForFigure(whichFigure)
% str = GetDataDirForFigure(whichFigure)
%
% Returns a string containing the data directory for the specified
% figure number. WHICHFIGURE must be in the range 2-8.

if (whichFigure<2 | whichFigure>8)
  error('Figure must be in the range 2-8.');
end

str = fullfile(GetRootDir('datafigs'),['figure' num2str(whichFigure)]);
