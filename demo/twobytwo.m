function P = twobytwo(varargin)
% TWOBYTWO A nice 2x2 parameterized matrix example
%
% Construct the 2x2 PME
%   [1+epsilon  t][x] = [2]
%   [t          1][y] = [1]
%
% Alternatively, construct the multi-dimensional 2x2 PME
%   [1+epsilon  t1] * ... * [1+epsilon  td][x] = [2]
%   [t1          1] * ... * [td          1][y] = [1]
%
% P = twobytwo() construct a problem where epsilon = 0.2
% P = twobytwo(epsilon) construct a problem with epsilon 
% P = twobytwo('epsilon',e) specify epsilon as a named argument
% P = twobytwo('dim',d) specify the dimensional as a named argument
%       dim=1 unless otherwise specified.
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
dim = [];
if nargin==1,
    epsilon = varargin{1};
else
    args = struct(varargin{:});
    if isfield(args,'epsilon'), epsilon=args.epsilon; args.epsilon = []; end
    if isfield(args,'dim'), dim=args.dim; args.dim = []; end
end
if isempty(epsilon), epsilon=0.2; end
if isempty(dim), dim=1; end

A0 = @(t) [1+epsilon t; t 1];
d = dim;

    function A=multidim_mat(s)
        if length(s)~=d, error('pmpack:wrongParameterSize',...
                'this 2x2 problem takes %i parameters, not %i', d, length(s)); end
        A = A0(s(1));
        for pi=2:d
            A = A*A0(s(pi));
        end
    end
        
P.b = @(t) [2; 1];
P.N = 2;        
P.d = 1;

% Create the anonymous functions for mat-vec interface.
if d==1
    P.A=@(t) [1+epsilon t; t 1];
    P.Av=@(t,v) P.A(t)*v;
    P.solve = @(t) [1+epsilon t; t 1]\[2;1];
    P.s = jacobi_parameter;
else
    P.A = @multidim_mat;
    P.Av = @(t,v) P.A(t)*v;    
    P.solve = @(t) P.A(t)\[2;1];
    variables(d) = jacobi_parameter;
    for i=1:d-1
        variables(i) = jacobi_parameter;
    end
    P.s = variables;
end

% main end
end
