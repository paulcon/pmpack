function P = twobytwo(varargin)
% TWOBYTWO A nice 2x2 parameterized matrix example
%
% Construct the 2x2 PME
%   [1+epsilon  t][x] = [2]
%   [t          1][y] = [1]
% The solution is
%   x = 
%
% P = twobytwo() construct a problem where epsilon = 0.2
% P = twobytwo(epsilon) construct a problem with epsilon 
% P = twobytwo('epsilon',e) specify epsilon as a named argument
%
% Return:
% The output P is a struct with the following fields:
% P.A - a parameterized matrix interface
% P.Av - a parameterized matrix-vector interface
% P.b - a parameterized right hand side
% P.N - the size of the matrix
% P.d - the number of parameters (1) for this problem
% P.solve - a parameterized solution
% P.s - the vector of parameters for this problem
% 
% Example:
%   P = twobytwo(0.01); % a harder problem
%   X = pseudospectral(P.solve,P.s,10);
%   residual_error_estimate(X,P.Av,P.b)
%   % eek, bad accuracy, use 100 points!
%   X = pseudospectral(P.solve,P.s,100);
%   residual_error_estimate(X,P.Av,P.b)

% Copyright, Stanford University, 2008-2010
% Paul G. Constantine, David F. Gleich

epsilon = [];
if nargin==1,
    epsilon = varargin{1};
else
    args = struct(varargin{:});
    if isfield(args,'epsilon'), epsilon=args.epsilon; args.epsilon = []; end
end
if isempty(epsilon), epsilon=0.2; end

% Create the anonymous functions for mat-vec interface.
P.A=@(t) [1+epsilon t; t 1];
P.Av=@(t,v) P.A(t)*v;
P.b=@(t) [2; 1];
P.N=2;
P.d = 1;
P.solve = @(t) [1+epsilon t; t 1]\[2;1];
P.s = jacobi_parameter;
