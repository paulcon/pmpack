function z=galerkin_preconditioner(R,x)
%GALERKIN_PRECONDITIONER
%
% z = galerkin_preconditioner(R,x)
%
% Used within the function handle to precondition the Galerkin system with
% the constant NxN matrix R. See examples for details.

N=size(R,1); n=length(x)/N;
X=reshape(x,N,n);
Z=R\X;
z=Z(:);
end