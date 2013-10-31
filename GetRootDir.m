% function rootDir = GetRootDir(whichRoot)
% 
% Returns the root directory specified. WHICHROOT is not case
% sensitive, and can be one of
%
% INSTALL: The top level directory in which /code and /data are
% contained.
%
% CODE: The top level code directory.
% 
% FIGURES: The top level code directory for the figures.
%
% DATAROOT: The root directory for the data, containing spike times,
% odors, and post-processed data.
%
% DATAPROC: The processed data directory, containing the results of
% computation on the spike times data.
%
% DATAFIGS: The subfolder in the processed data directory containing
% the results of computations by the PROCESSDATA scripts for each
% figure.
% 
% SPIKETIMES: The directory containing the PN and KC spike times.
%
% ODORS: The directory containing information about the odors used.

function rootDir = GetRootDir(whichRoot)
installRoot = INSTALL_ROOT;

if (nargin==0)
  whichRoot = 'install';
end

switch (lower(whichRoot))
 case 'install'
  rootDir = installRoot;
 case 'code'
  rootDir = fullfile(installRoot, 'code');
 case 'figures'
  rootDir = fullfile(installRoot, 'code', 'figures');
 case 'dataroot'
  rootDir = fullfile(installRoot, 'data', 'figures');
 case 'dataproc'
  rootDir = fullfile(installRoot, 'data', 'proc');
 case 'datafigs'
  rootDir = fullfile(installRoot, 'data', 'proc', 'figures');
 case 'spiketimes'
  rootDir = fullfile(installRoot, 'data', 'spt');
 case 'odors'
  rootDir = fullfile(installRoot, 'data', 'odors');
 otherwise
  error('Unknown root directory for "%s".', whichRoot);
end
