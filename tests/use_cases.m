% use case script

%% clean up
close all; clear all;
if matlabpool('size'), matlabpool('close'); end

%% define 2x2 test problem with one parameter
epsilon=0.5;
A=@(s) [1+epsilon s; s 1];
Ax=@(s,x) A(s)*x;
b=@(s) [2; 1];
iAb=@(s) A(s)\b(s);
N=2;
s=parameter();

%% Standard pseudospectral
p_order=15;
X=pseudospectral(iAb,s,p_order);

%% pseudospectral convergence studies
[X,errz]=pseudospectral(iAb,s,'adapt'); % using default relative error and tolerance
[X,errz]=pseudospectral(iAb,s,'adapt','ptol',1e-6); % high tolerance
convtype='mincoeff';
[X,errz]=pseudospectral(iAb,s,'adapt','ConvType',convtype);

%% convergence study with residual error estimate
convtype='resid';
[X,errz]=pseudospectral(iAb,s,'adapt','ConvType',convtype,'matfun',A,'vecfun',b);

%% convergence study with a reference solution
Y=pseudospectral(iAb,s,50);
[X,errz]=pseudospectral(iAb,s,'adapt','RefSoln',Y);

%% Standard Galerkin
% Builds full Galerkin matrix and solves with backslash.
p_order=15;
X=spectral_galerkin(A,b,s,p_order);

%% Standard Galerkin with high order integration
X=spectral_galerkin(A,b,s,p_order,'qorder',4*p_order);

%% Low mem option 
% like an unassembled finite element method
X=spectral_galerkin(A,b,s,p_order,'Solver',@minres,'lowmem',1);

%% Matvec option with optional inputs for pcg
pcgtol=1e-10; pcgmaxi=1000;
solver=@(A,b) pcg(A,b,pcgtol,pcgmaxi);
X=spectral_galerkin(Ax,b,s,p_order,'Solver',solver);

%% Preconditioning
s_mid=midpoint(s);
R=chol(A(s_mid),'lower');
pcon1=@(x) galerkin_preconditioner(R',x);
pcon2=@(x) galerkin_preconditioner(R,x);
solver=@(A,b) pcg(A,b,pcgtol,pcgmaxi,pcon1,pcon2);
X=spectral_galerkin(Ax,b,s,p_order,'Solver',solver);

%% Convergence study with preconditioning
[X,errz]=spectral_galerkin(Ax,b,s,'adapt','Solver',solver,'ConvType','resid');