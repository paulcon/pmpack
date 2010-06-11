% galerkin test script

%% Setup

%%
% Problem 1
iAb=@(p) p^4;
A=@(p) 1; 
b=@(p) p^4;
Ax=@(p,x) x;
s=parameter();

%%
% Problem 2
iAb2=@(p) p(1)^4*p(2)^9;
A2=@(p) 1; 
b2=@(p) p(1)^4*p(2)^9+2;
Ax2=@(p,x) x;
s2=[parameter(); parameter()];

%% Compute exact solutions
Y=pseudospectral(iAb,s,10);
Y2=pseudospectral(iAb2,s2,10);

%% Test basic galerkin
[X,errz]=spectral_galerkin(A,b,s,5);
[X,errz]=spectral_galerkin(A2,b2,s2,[5 5]);
[X,errz]=spectral_galerkin(Ax2,b2,s2,[5 7]);
[X,errz]=spectral_galerkin(A2,b2,s2,[5 7],'qorder',[4 4]);
%% Test adaptive methods
[X,errz]=spectral_galerkin(A,b,s,'adapt');
[X,errz]=spectral_galerkin(A,b,s,'adapt','ErrEst','MinCoeff');
[X,errz]=spectral_galerkin(A2,b2,s2,'adapt','ErrEst','MinCoeff');
[X,errz]=spectral_galerkin(A,b,s,'adapt','ErrEst','RelErr');
[X,errz]=spectral_galerkin(A2,b2,s2,'adapt','ErrEst','RelErr');
%% Test reference solutions
[X,errz]=spectral_galerkin(A,b,s,'adapt','RefSoln',Y);
[X,errz]=spectral_galerkin(A2,b2,s2,'adapt','RefSoln',Y2,'ptol',1e-13);
%%
[X,errz]=spectral_galerkin(Ax2,b2,s2,'adapt','solver',@pcg);
[X,errz]=spectral_galerkin(A2,b2,s2,'adapt','solver',@pcg,'lowmem',1);
[X,errz]=spectral_galerkin(Ax2,b2,s2,'adapt','solver',@pcg,'lowmem',1);
[X,errz]=spectral_galerkin(A2,b2,s2,'adapt','solver',@pcg,'lowmem',1,'ptol',1e-13);
