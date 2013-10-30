function [Y, odorNames] = GetOdorNamesAsBinaryVectors(varargin)
% [Y, odorNames] = GetOdorNamesAsBinaryVectors(varargin)
%
% Returns the odors list as a 40 x 8 matrix whose rows are the odors as binary vectors, in ABCDWXYZ format.
% The *high odors and paraffin oil are not returned. The list of odors included is returned in ODORNAMES.
% An optional argument of 'full' will return the full odors list as binary vectors.
load(fullfile(GetRootDir('odors'),'odorschar.mat'));
load(fullfile(GetRootDir('odors'),'odors.mat'));

if (~isempty(varargin) & isequal(lower(varargin{1}),'full'))
  whichOdors = [1:44];
else
  whichOdors = [4:11 13:44];
end

X = odorschar;
X = X(whichOdors,:);
cmps = 'ABCDWXYZ';
Y = 0*X;
for i = 1:8
  Y(:,i) = arrayfun(@(j) sum(X(j,:)==double(cmps(i)))>0,1:numel(whichOdors))';
end
odorNames = odors(whichOdors);

