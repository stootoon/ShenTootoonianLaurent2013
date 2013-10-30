function name = odor_index_to_name(odorIndex, varargin)
% name = odor_index_to_name(odorIndex, varargin)
%
% Returns the odor name (in the complex mixture experiments) with index
% ODORINDEX.  An optional second argument specifies a file from which
% contains the odor names. Each must be prefixed by 'odor'.
%
% Example: Return the odor with index 24:
%
% >> odor_index_to_name(24)
%
% ans =
% 
% odorYZ

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


if (nargin>1) % An odors list file was provided, try to use it.
  odorFile = varargin{1};
  
  if (~fileexists(odorFile))
    error('Could not open odor file "%s" for read.',odorFile);
  end
  
  odorOrder = {};
  nOdors = 1;
  fp = fopen(fileExists,'r');
  
  while(~feof(fp)) % Read strings from the file and store those prefixed by 'odor'
    odor = fscanf('%s',fp);
    if (strmatch('odor', odor)==1)
      odorOrder{nOdors} = odor;
      nOdors = nOdors+1;
    end
  end
  
  disp(['Read ' num2str(nOdors) ' from ' odor '.']);
  fclose(fp);
end

if (odorIndex<1 || odorIndex>length(odorOrder))
    error('Odor index (%d) is out of range.');
end

name = odorOrder{odorIndex};
   
