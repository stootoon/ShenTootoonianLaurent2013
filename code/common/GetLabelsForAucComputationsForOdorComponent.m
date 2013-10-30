function L = GetLabelsForAucComputationsForOdorComponent(whichCmp)
% L = GetLabelsForAucComputationsForOdorComponent(whichCmp)
%
% Returns a vector L such that L(i) = 1 if odors{i} contains the 
% specified component, and -1 otherwise, and 0 for paraffin oil.
%
% Example:
%
% L = GetIngroupLabelsForOdorComponent('A');
%
% L(odor_name_to_index('odorA'))
% L(4)
%
% ans =
%
%     1

odors = GetOdorsList;
odors{44} = 'odorABCDWXYZ';
L = cellfun(@(x) ~isempty(strfind(x(4:end), whichCmp)), odors);
L = double(L);
L = 2*L-1;
L(12) = 0; % Paraffin oil is always out of the group.


