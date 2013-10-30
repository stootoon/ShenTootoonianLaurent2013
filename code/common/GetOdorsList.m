function odors = GetOdorsList()
% odors = GetOdorsList()
%
% Returns the list of odors used in Kai's experiments.

odors = LoadVarFromMatFileByName(fullfile(GetRootDir('odors'),'odors.mat'),'odors');




