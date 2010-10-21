function P = twobytwo_func(epsilon)
% TWOBYTWO A 2x2 parameterized matrix example
%
% Construct the 2x2 PME
%   [1+epsilon  s][x] = [2]
%   [s          1][y] = [1]
%
% Alternatively, construct the multi-dimensional 2x2 PME
%   [1+epsilon1  t1] * ... * [1+epsilond  td][x] = [2]
%   [t1           1] * ... * [td           1][y] = [1]
%
% P = twobytwo(epsilon) construct a problem with the d-dimensional array
% epsilon.
%
% Return:
% The output P is a struct with the following fields:
% P.A - a parameterized matrix interface
% P.Av - a parameterized matrix-vector interface
% P.b - a parameterized right hand side
% P.N - the size of the matrix
% P.d - the number of parameters for this problem
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

d = length(epsilon);

A0 = @(t,e) [1+e t; t 1];

function A=multidim_mat(s)
    if length(s)~=d, error('Length of s=%i but length of epsilon=%i',length(s),d); end
    A = A0(s(1),epsilon(1));
    for pi=2:length(s)
        A = A*A0(s(pi),epsilon(pi));
    end
end
        
P.b = @(t) [2; 1];
P.N = 2;        
P.d = d;

% Create the anonymous functions for mat-vec interface.
P.A = @multidim_mat;
P.Av = @(t,v) P.A(t)*v;    
P.solve = @(t) P.A(t)\[2;1];
variables(d) = jacobi_parameter;
for i=1:d-1
    variables(i) = jacobi_parameter;
end
P.s = variables;

% main end
end
