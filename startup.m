function startup()
% startup()
%
% Startup script for ShenTootoonianLaurent2013.
% Adds the core code directories to the path.

installRoot = GetRootDir;
addpath(installRoot);
addpath(fullfile(installRoot, 'code'));
addpath(fullfile(installRoot, 'code', 'common'));
addpath(genpath(fullfile(installRoot, 'code', 'util')));
