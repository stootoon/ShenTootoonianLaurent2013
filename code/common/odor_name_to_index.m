function idx = odor_name_to_index(odorName, varargin)
% idx = odor_name_to_index(odorName, varargin)
%
% Returns the index of the odor ODORNAME (e.g. 'odorYZ') in the list
% of odors in the complex mixtures experiments, or -1 if the odor
% could not be found. An optional second argument can provide the name
% of a file containing an alternate list of odors to use. Each must be
% prefixed by 'odor'.
%
% Example:
%
% odor_name_to_index('odorYZ')
% 
% ans =
%
%    24
%
% odor_name_to_index('odorFGH')
% 
% ans =
%
%    -1

odorOrder = {
'odorAhigh',
'odorChigh',
'odorWhigh',
'odorA',
'odorB',
'odorC',
'odorD',
'odorW',
'odorX',
'odorY',
'odorZ',
'odorP',
'odorAB',
'odorAC',
'odorAD',
'odorAX',
'odorBC',
'odorDW',
'odorWX',
'odorWY',
'odorWZ',
'odorXY',
'odorXZ',
'odorYZ',
'odorABC',
'odorACD',
'odorAXZ',
'odorBCX',
'odorBDW',
'odorDXY',
'odorWXY',
'odorWYZ',
'odorABCD',
'odorABCX',
'odorBCWX',
'odorBDWX',
'odorDWYZ',
'odorWXYZ',
'odorABCDX',
'odorABCWX',
'odorADWYZ',
'odorAWXYZ',
'odorBCWXZ',
'odorALL'};

if (nargin>1) % An alternate odor list file was provided, try to use it.
  
  odorFile = varargin{1};
  
  if (~fileexists(odorFile))
    error('Could not open odor file "%s" for read.',odorFile);
  end
  
  odorOrder = {};
  nOdors = 1;
  fp = fopen(fileExists,'r');
  
  while(~feof(fp)) % Read the odor names, and store those prefixed by 'odor'.
    odor = fscanf('%s',fp);
    if (strmatch('odor',odor)==1)
      odorOrder{nOdors} = odor;
      nOdors = nOdors+1;
    end
  end
  
  disp(['Read ' num2str(nOdors) ' from ' odor '.']);
  fclose(fp);

end

% Now look for the specified odor in the list and return its index.
idx = -1;
for i=1:length(odorOrder)
  if (strcmp(odorName, odorOrder{i}))
    idx = i;
    break;
  end
end

