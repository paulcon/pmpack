function P=pmpack_problem(pname,varargin)
% PMPACK_PROBLEM Load a problem structure for a parameterized matrix problem
%
% P = pmpack_problem(pname);
%
% Outputs:
%   P: A problem struct where 
%           P.A : is a function returning the parameterized matrix
%           P.b : is a function returning the parameterized rhs
%           P.s : is the set of parameters
%           P.N : is the size of the matrix
%           P.solve : is a function to solve P.A(s0)\P.b(s0) at a point s0
% Inputs:
%   pname: the name of a problem in the problem directory
%
% As a side effect, this function changes the Matlab path to include the
% problems directory of pmpack.  
%
% Example:
%   P = pmpack_problem('twobytwo');
%
% See also TWOBYTWO NOBILE_TEMPONE_WEBSTER ELLIPTIC1D

% Copyright 2009-2010 David F. Gleich (dfgleic@sandia.gov) and Paul G. 
% Constantine (pconsta@sandia.gov)
%
% History
% -------
% :2010-06-14: Initial release


fullpath = mfilename('fullpath');
filepath = fileparts(fullpath);
probdir = fullfile(filepath,'problems');
addpath(probdir);
P = feval(pname,varargin{:});
