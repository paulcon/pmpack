function z=galerkin_preconditioner(R,x)
N=size(R,1); n=length(x)/N;
X=reshape(x,N,n);
Z=R\X;
z=Z(:);
end