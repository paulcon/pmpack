function P = elliptic1d(varargin)
% ELLIPTIC1D A simple elliptic ODE as a parameterized matrix problem
%
% Consider the following ODE:
%       d/dt (a(s,t)*du/dt) = 1, t\in [0,1]
%       u(0)=u(1)=0
% where 
%       a(s,t)=1+4*s*(t^2-t), s\in [-1,1-epsilon]
% This function sets up the matrix calls and matrix-vector multiply interface
% for the discretized version of this problem using a finite element method
% with a standard piece-wise linear basis. 
%
% P = elliptic1d() construct a problem where the matrix has size 511 and
%   epsilon=0.01
% P = elliptic1d(epsilon) construct a problem where the matrix has size
%   511 and epsilon is given by the parameter
% P = elliptic1d('size',N,'epsilon',e) specificy both parameters
%
% Return:
%
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
%   P = elliptic1d('epsilon',0.01); errs = zeros(19,1);
%   for n=2:20 % try pseudospectral with 2 to 10 basis functions
%     X = pseudospectral(P.solve,P.s,n);
%     errs(n-1) = size(X,1)*residual_error_estimate(X,P.Av,P.b)/norm(P.b(0));
%   end
%   semilogy(2:20,errs,'.-');
% 
% See also ELLIPTICND

% Copyright, Stanford University, 2009-2010
% Paul G. Constantine, David F. Gleich

Nel = [];
epsilon = [];
if nargin==1,
    epsilon = varargin{1};
else
    args = struct(varargin{:});
    if isfield(args,'size'), Nel=args.size; args.size = []; end
    if isfield(args,'epsilon'), epsilon=args.epsilon; args.epsilon = []; end
end
if isempty(Nel), Nel=512; end
if isempty(epsilon), epsilon=0.01; end

% Set up the 1d grid.
grid.Nel=Nel;
grid.h=1/grid.Nel;
grid.Ndof=grid.Nel-1;
grid.x=0:grid.h:1;

% Integrated form of the elliptic coefficients for the finite element
% stiffness matrix.
A1=@(x) 4*((x.^3/3) - (x.^2/2));
A0=1;
f0=1;

% Construct the stiffness matrices and the right hand side.
K0=A0*grid.h*full(gallery('tridiag',grid.Ndof,-1,2,-1));
K1=stiffness(A1,grid);
f=f0*grid.h^3*ones(grid.Ndof,1);

% Create the anonymous functions for mat-vec interface.
% Notice that the parameteric input has been scaled to take the stanadard
% interval [-1,1].
A=@(t) sparse(K0)+(0.5*(t+1)*(2-epsilon)-1)*sparse(K1);
Av=@(t,v) sparse(K0)*v+(0.5*(t+1)*(2-epsilon)-1)*sparse(K1)*v;
b=@(t) f;

% initialize the size of th matrix.
N=grid.Ndof;

P.N = size(K0,1);
P.A = A;
P.Av = Av;
P.b = b;
P.d = 1;
P.b_const = 1; 
P.solve = @(t) A(t)\f;
P.s = jacobi_parameter();

function r = stiffness(A,grid)

Nel=grid.Nel;
Ndof=grid.Ndof;
x=grid.x;

K=zeros(Ndof);
K(1,1)=A(x(3))-A(x(1));
K(1,2)=-A(x(3))+A(x(2));
for i=2:Ndof-1
    K(i,i)=A(x(i+2))-A(x(i));
    K(i,i+1)=-A(x(i+2))+A(x(i+1));
    K(i,i-1)=-A(x(i+1))+A(x(i));
end
K(Ndof,Ndof)=A(x(Nel+1))-A(x(Nel-1));
K(Ndof,Ndof-1)=-A(x(Nel))+A(x(Nel-1));

r = K;

