function parent = GetCurrentDirFromPathString(str)
% parent = GetCurrentDirFromPathString(str)
% 
% Returns the name of the current directory from the path string STR
% as the text following the last file separator character.

indLastSep = max(find(str==filesep));

if (isempty(indLastSep))
  parent = str;
elseif (indLastSep == numel(str))
  parent = '';
else
  parent = str(indLastSep+1:end);
end

