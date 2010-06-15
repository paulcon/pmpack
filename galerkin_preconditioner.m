function z=galerkin_preconditioner(R,x)
%GALERKIN_PRECONDITIONER Apply one matrix to precondition a Galerkin system
%
% z = galerkin_preconditioner(R,x)
%
% Compute the action z = kron(speye(n),R^{-1})*x in an efficient manner.
% The size of x must be size(R)*n.  This function is used to precondition
% the Galerkin system with a single NxN matrix R.  See the example for its
% use.
%
% Example:  
%   P = pmpack_problem('elliptic1d');
%   s0 = midpoint(P.s);
%   R = chol(P.A(s0));
%   P1 = @(x) galerkin_preconditioner(R',x);
%   P2 = @(x) galerkin_preconditioner(R,x);
%   X1 = spectral_galerkin(P.A,P.b,P.s,2,...
%         'solver',@(A,b) minres(A,b,1e-5,100,P1,P2));
%   X2 = spectral_galerkin(P.A,P.b,P.s,2,... % no preconditioning
%         'solver',@(A,b) minres(A,b,1e-5,1500));
%
% See also SPECTRAL_GALERKIN

% Copyright 2009-2010 David F. Gleich (dfgleic@sandia.gov) and Paul G. 
% Constantine (pconsta@sandia.gov)
%
% History
% -------
% :2010-06-14: Initial release

N=size(R,1); n=length(x)/N;
X=reshape(x,N,n);
Z=R\X;
z=Z(:);
