function varargout = MatlabPoolWrapper(varargin)
% varargout = MatlabPoolWrapper(varargin)
%
% A wrapper for MATLABPOOL. Checks to see if MATLABPOOL exists, and
% if so, will run it transparently. Otherwise, will do nothing.
%
% We only every use MATLABPOOL to open, close, check size, so it's ok
% if we let these fall through.

try 
  matlabpool('size');
catch % matlabpool doesn't work, so just let 
  varargout{1:nargout} = nan*(1:nargout);
  return;
end
[varargout{1:nargout}] = matlabpool(varargin{:});
