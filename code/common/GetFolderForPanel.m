function folder = GetFolderForPanel(whichPanel)
% folder = GetFolderForPanel(whichPanel)
%
% Returns the folder containing the ProcessData script producing data
% for the specified panel, unless the panel is generated directly, in
% which case the filename of the generating file is returned.

switch upper(whichPanel)
 case {'2A','2B','2C','2D'}
  folder = fullfile(GetCodeDirForFigure(2), 'MakeBinaryMixtureRasterExample.m');
 case {'2E', '2F', 'S2A', 'S2B', 'S2C', 'S2D', 'S2E', 'S2F'}
  folder = fullfile(GetCodeDirForFigure(2), 'BinaryMixtureResponseLinearity');
 case '3A'
  disp('LLE not yet implemented.');
  folder = '';
 case {'3B', '3C'}
  folder = fullfile(GetCodeDirForFigure(3), 'CorrelationDistanceMatrices');
 case {'3D', '3E', '3F'}
  disp('Correlation matrices available. LLE not yet implemented.');
  folder = fullfile(GetCodeDirForFigure(3), 'CorrelationDistanceMatrices');
 case {'3G', 'S3A', 'S3B'}
  folder = fullfile(GetCodeDirForFigure(3), 'BinaryMixtureTrajectoryProjections');
 case '3H'
  folder = fullfile(GetCodeDirForFigure(3), 'BinaryMixtureTrajectoryClustering');
 case {'3I', 'S3C', 'S3D'}
  folder = fullfile(GetCodeDirForFigure(3), 'BinaryMixtureTrajectoryMorphs');
 case {'4A', '4B', 'S4A', 'S4B', 'S4C', 'S4D', 'S4E', 'S4F', 'S4G'}
  folder = GetCodeDirForFigure(4);
 case {'5A', '5B', '5C'}
  disp('Correlation matrices available. LLE not yet implemented.');
  folder = fullfile(GetCodeDirForFigure(5), 'CorrelationDistanceMatrices');
 case {'5D', 'S5A', 'S5B'}
  folder = fullfile(GetCodeDirForFigure(5), 'ComplexMixtureGlobalTrajectoryDistanceWithOverlap_noduplicates');
 case '5E'
  folder = fullfile(GetCodeDirForFigure(5), 'ComplexMixtureTransientProximity_noduplicates');
 case {'5F', 'S5C', 'S5D'}
  folder = fullfile(GetCodeDirForFigure(5), 'ComplexMixturePerBinTrajectoryDistanceWithOverlap_noduplicates');
 case {'6A', '6B', '6C', '6D'}
  folder = fullfile(GetCodeDirForFigure(6), 'MakeFigures.m');
 case {'6E', '6F'}
  folder = fullfile(GetCodeDirForFigure(6), 'PnKcAuc_noduplicates');
 case {'S6A', 'S6B', 'S6C', 'S6D'}
  folder = fullfile(GetCodeDirForFigure(6), 'PnAndKcPromiscuity_noduplicates','MakeFiguresForPaper.m');
 case {'7A'}
  folder = fullfile(GetCodeDirForFigure(7), 'KcResponseLatency');
 case {'7B','7C','7D', '7E'}
  fprintf('Panel "%s" is not yet available.\n', upper(whichPanel));
 case {'7F','7I'}
  fprintf('Panel "%s" is a schematic.\n', upper(whichPanel));
 case {'7G', '7H'}
  folder = fullfile(GetCodeDirForFigure(7), 'PnKcConnectivity2_noduplicates');
 case {'7J', '7K'}
  folder = fullfile(GetCodeDirForFigure(7), 'ComplexMixtureTrajectoryReconstruction');
 case 'S7A'
  folder = fullfile(GetCodeDirForFigure(7), 'PnKcResponseTimecourse_noduplicates');
 case {'8A', '8B', '8C', 'S8A'}
  folder = fullfile(GetCodeDirForFigure(8), 'classification_noduplicates');
 case {'S8B', 'S8C'}
  folder = fullfile(GetCodeDirForFigure(8), 'generalization_noduplicates');
 case {'8D','8E','8F','8G','8I'}
  fprintf('Panel "%s" is a schematic.\n', upper(whichPanel));
 otherwise
  fprintf('Panel "%s" is not available.\n', upper(whichPanel));
end
  