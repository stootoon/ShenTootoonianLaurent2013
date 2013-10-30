function g = OdorGroupingFunction1(odorList)
% g = OdorGroupingFunction1(odorList)
% 
% Assigns group ids to odor names. Odors must be prefixed by 'odor'.
% -1: Paraffin Oil
% 0: 'high' odors
% 8: ALL odors
% 1:7: Assigned according to number of components

g = nan(1,numel(odorList));

for i = 1:numel(odorList)
    od = odorList{i};
    if (od(end-3:end)=='high')
        g(i)=0;
    elseif (od(end)=='P')
        g(i) = -1;
    elseif (od(end-2:end)=='ALL')
        g(i) = 8;
    else
        g(i) = length(od)-4;
    end
end