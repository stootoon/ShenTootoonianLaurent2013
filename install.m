function success = install()
% function success = install()
% 
% Installs the code for ShenTootoonianLaurent2013 in the current
% directory by 
%
% 1) Setting INSTALL_ROOT in code/common/GetRootDir.m to the current
% directory, and
%
% 2) Compiling code/util/moddecomp.c 

success = false;
%% Update INSTALL_ROOT in GetRootDir.m
% Try to open 'GetRootDir'
getRootDirFile = './GetRootDir.m';
[fp, errMsg] = fopen(getRootDirFile, 'r');
if (fp == -1)
  fprintf('Could not open "%s" for read: %s\n', getRootDirFile, errMsg);
  fprintf('Installation aborted.\n\n');
  return;
end
fclose(fp);

% File exists, so look for the tag to replace...
contents = fileread(getRootDirFile);
indRep   = strfind(contents, 'INSTALL_ROOT');
if (isempty(indRep))
  fprintf('Could not find "INSTALL_ROOT" in "%s".\n', getRootDirFile);
  fprintf('Installation aborted.\n\n');
  return;
end

% INSTALL_ROOT was found, so replace it.
installRoot = pwd;
contents = strrep(contents, 'INSTALL_ROOT', ['''' installRoot '''']);

% Write the new contents to file.
[fp, errMsg] = fopen(getRootDirFile, 'w');
if (fp == -1)
  fprintf('Could not open "%s" for write: %s\n', errMsg);
  fprintf('Installation aborted.\n\n');
  return;
end
fprintf(fp, '%s', contents);
fclose(fp);

fprintf('Set INSTALL_ROOT to "%s".\n', installRoot);

%% Compute code/util/moddecomp.c
retVal = mex('code/util/moddecomp.c', '-outdir', 'code/util');
if (retVal)
  fprintf('Could not compile "code/util/moddecomp.c"\n');
  fprintf('Installation aborted.\n\n');
  return;
end
fprintf('Compiled "code/util/moddecomp.c".\n');

%% All done.
fprintf('\nInstallation successful.\n\n');
success = true;
return;
