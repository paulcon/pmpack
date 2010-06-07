function P=ellipticnd(varargin)
% ELLIPTICND A simple elliptic ODE as a 5d parameterized matrix problem
%
% Consider the following ODE:
%       d/dt (a(s,t)*du/dt) = 1, t\in [0,1]
%       u(0)=u(1)=0
% where a(s,t) is a log normal approximation as given in Nobile, et al.
% This script sets up the matrix calls and matrix-vector multiply interface
% for the discretized version of this problem using a finite element method
% with a standard piece-wise linear basis. 
%
% P = ellipticnd() construct a problem where the matrix has size 511 
% P = ellipticnd('size',N) construct a problem where the matrix has size N
%
% Return:
%
% The output P is a struct with the following fields:
% P.A - a parameterized matrix interface
% P.Av - a parameterized matrix-vector interface
% P.b - a parameterized right hand side
% P.N - the size of the matrix
% P.d - the number of parameters (5) for this problem
% P.solve - a parameterized solution
% P.s - the vector of parameters for this problem
% 
% Example:
%   P = ellipticnd;
%   dt=tic; X1 = pseudospectral(P.solve,P.s,1); toc(dt);
%   dt=tic; X2 = pseudospectral(P.solve,P.s,2); toc(dt); % a lot longer!
% 
% See also ELLIPTIC1D

% Copyright, Stanford University, 2009-2010
% Paul G. Constantine, David F. Gleich

Nel = [];
if nargin==1,
    epsilon = varargin{1};
else
    args = struct(varargin{:});
    if isfield(args,'size'), Nel=args.size; args.size = []; end
end
if isempty(Nel), Nel=512; end


% Set number of elements
grid.Nel=Nel;
grid.h=1/grid.Nel;
grid.Ndof=grid.Nel-1;
grid.x=0:grid.h:1;

% Create the anonymous functions for mat-vec interface.
A=@(t) stiffness_system(t,grid);
Av=@(t,v) A(t)*v;

% Create anonymous function for right hand side.
f0=1;
f=f0*grid.h^2*cos(grid.x(2:end-1)); f=f(:);
b=@(t) f;

% initialize the size of th matrix.
N=grid.Ndof;

P.N = N;
P.A = A;
P.Av = Av;
P.b = b;
P.d = 5;
P.b_const = 1; 
P.solve = @(t) A(t)\f;
P.s = [jacobi_parameter() jacobi_parameter() jacobi_parameter() jacobi_parameter() jacobi_parameter()];


function K = stiffness_system(s,grid)

A0=0.5; 
K0=A0*gallery('tridiag',grid.Ndof,-1,2,-1);

Lc=0.5; Lp=max(1,2*Lc); L=Lc/Lp;
spi=sqrt(sqrt(pi)*L/3);
A1=@(x)exp(1 + ...
    (s(1)*spi/sqrt(2)) + ...
    spi*exp(-(pi*L)^2/8)*sin(pi*x/Lp)*s(2) + ...
    spi*exp(-(pi*L)^2/8)*cos(pi*x/Lp)*s(3) + ...
    spi*exp(-(2*pi*L)^2/8)*sin(2*pi*x/Lp)*s(4) + ...
    spi*exp(-(2*pi*L)^2/8)*cos(2*pi*x/Lp)*s(5));

K1=sparse(grid.Ndof,grid.Ndof);
K1(1,1)=2*A1(grid.x(2));
K1(1,2)=-A1(0.5*(grid.x(3)+grid.x(2)));
for i=2:grid.Ndof-1
    K1(i,i)=2*A1(grid.x(i+1));
    K1(i,i+1)=-A1(0.5*(grid.x(i+2)+grid.x(i+1)));
    K1(i,i-1)=-A1(0.5*(grid.x(i+1)+grid.x(i)));
end
K1(grid.Ndof,grid.Ndof)=2*A1(grid.x(grid.Nel));
K1(grid.Ndof,grid.Ndof-1)=-A1(0.5*(grid.x(grid.Nel)+grid.x(grid.Nel-1)));

K=K0+K1;
