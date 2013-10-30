function L = GetBalancedIngroupLabelsForOdorComponent(whichCmp)
% L = GetBalancedIngroupLabelsForOdorComponent(whichCmp)
%
% Returns a vector L such that L(i) = 1 if odors{i} contains the
% specified component, -1 if not, and 0 if not included. The labeling
% is balanced following Kai's description in the methods so that every
% outgroup element corresponds to an ingroup element but with the
% component in question removed. Hence some odors are not used and
% will be assigned the label 0.
%
% Example:
%
% L = GetBalancedIngroupLabelsForOdorComponent('A');
%
% L(odor_name_to_index('odorA'))
% L(4)
%
% ans =
%
%     1

load(fullfile(GetRootDir('odors'),'odors.mat'));
L = zeros(1,44);
switch(whichCmp)
 case 'A'
  inOdors = 'Ahigh,A,AB,AC,AD,AX,ABC,ACD,AXZ,ABCD,ABCX,ABCDX,ABCWX,ADWYZ,AWXYZ';
  outOdors= 'Chigh,B,C,D,X,BC,DW,XZ,BCX,BDW,DXY,WXY,WYZ,BCWX,BDWX,DWYZ,WXYZ,BCWXZ';
 case 'B'
  inOdors = 'B,AB,BC,ABC,BCX,BDW,ABCD,ABCX,BCWX,BDWX,ABCDX,ABCWX,BCWXZ';
  outOdors= 'Ahigh,A,Chigh,C,X,AC,AD,AX,DW,WX,ACD,AXZ,DXY,WXY,DWYZ,WXYZ,ADWYZ,AWXYZ';
 case 'C'
  inOdors = 'Chigh,C,AC,BC,ABC,ACD,BCX,ABCD,ABCX,BCWX,ABCDX,ABCWX,BCWXZ';
  outOdors= 'Ahigh,Whigh,A,B,X,AB,AD,WX,AXZ,BDW,DXY,BDWX,DWYZ,WXYZ,ADWYZ,AWXYZ';
 case 'D'
  inOdors = 'D,AD,DW,ACD,BDW,DXY,ABCD,BDWX,DWYZ,ABCDX,ADWYZ';
  outOdors= 'Whigh,A,B,W,AC,XY,ABC,WXY,WYZ,ABCX,BCWX,WXYZ,ABCWX,AWXYZ,BCWXZ';
 case 'W'
  inOdors = 'Whigh,W,DW,WX,WY,WZ,BDW,WXY,WYZ,BCWX,BDWX,DWYZ,WXYZ,ABCWX,BCWXZ';
  outOdors= 'Ahigh,Chigh,B,D,X,Y,Z,BC,XY,XZ,YZ,ABC,ACD,AXZ,BCX,DXY,ABCD,ABCX,ABCDX';
 case 'X'
  inOdors = 'X,AX,WX,XY,XZ,AXZ,BCX,DXY,WXY,ABCX,BCWX,BDWX,WXYZ,ABCDX,AWXYZ';
  outOdors = 'Ahigh,Chigh,Whigh,A,D,W,Y,Z,BC,WY,WZ,YZ,ABC,ACD,BDW,WYZ,ABCD,DWYZ,ADWYZ';
 case 'Y'
  inOdors = 'Y,WY,XY,YZ,DXY,WXY,WYZ,DWYZ,WXYZ,ADWYZ,AWXYZ';
  outOdors = 'D,W,X,Z,DW,WX,WZ,XZ,AXZ,BDW,BCWX,BDWX,ABCWX,BCWXZ';
 case 'Z'
  inOdors = 'Z,WZ,XZ,YZ,AXZ,WYZ,DWYZ,WXYZ,ADWYZ,AWXYZ,BCWXZ';
  outOdors = 'W,X,Y,AX,DW,XY,WY,BDW,WXY,DXY,BCWX,BDWX,ABCDX,ABCWX';
 otherwise 
  error('Expected the desired component to be one of {A, B, C, D, W, X, Y, Z}.');
end

inOdors = regexp(inOdors,',','split');
outOdors = regexp(outOdors,',','split');

inOdorInds = cellfun(@(x) odor_name_to_index(['odor' x]), inOdors);
outOdorInds= cellfun(@(x) odor_name_to_index(['odor' x]), outOdors);

L(inOdorInds) = 1;
L(outOdorInds) = -1;





