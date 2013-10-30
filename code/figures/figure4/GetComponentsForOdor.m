function indCmps = GetComponentsForOdor(whichOdor)
% indCmps = GetComponentsForOdor(whichOdor)
%
% Returns a vector of indices from 1-8 indicating the components
% presented in the specified odor. WHICHODOR can be a string, (either
% of the form 'odorAB', or 'AB'), or an integer from 1-44, indicating
% the odor id.

if (ischar(whichOdor))
  if (length(whichOdor)>4 && isequal(whichOdor(1:4),'odor'))
    odorId = odor_name_to_index(whichOdor);
  else
    odorId = odor_name_to_index(['odor' whichOdor]);
  end
  if (odorId < 1)
    error('Could not find odor "%s".', whichOdor);
  end
else
  odorId = whichOdor;
  if (odorId<1 || odorId>44)
    error('Odor ID out of range (must be in [1-44]).');
  end
end

B = GetOdorNamesAsBinaryVectors('full');
indCmps = find(B(odorId,:));


