% galerkin test script

iAb=@(p) p^4;
A=@(p) 1; 
b=@(p) p^4;
Ax=@(p,x) x;
s=parameter();
iAb2=@(p) p(1)^4*p(2)^9;
A2=@(p) 1; 
b2=@(p) p(1)^4*p(2)^9;
Ax2=@(p,x) x;
s2=[parameter(); parameter()];

Y=pseudospectral(iAb,s,10);
Y2=pseudospectral(iAb2,s2,10);

[X,errz]=spectral_galerkin(A,b,s,5);
[X,errz]=spectral_galerkin(A2,b2,s2,[5 5]);
[X,errz]=spectral_galerkin(Ax2,b2,s2,[5 7]);
[X,errz]=spectral_galerkin(A2,b2,s2,[5 7],'qorder',[4 4]);
[X,errz]=spectral_galerkin(A,b,s,'adapt');
[X,errz]=spectral_galerkin(A,b,s,'adapt','ConvType','MinCoeff');
[X,errz]=spectral_galerkin(A2,b2,s2,'adapt','ConvType','MinCoeff');
[X,errz]=spectral_galerkin(A,b,s,'adapt','ConvType','RelErr');
[X,errz]=spectral_galerkin(A2,b2,s2,'adapt','ConvType','RelErr');
[X,errz]=spectral_galerkin(A,b,s,'adapt','RefSoln',Y);
[X,errz]=spectral_galerkin(A2,b2,s2,'adapt','RefSoln',Y2,'ptol',1e-13);

[X,errz]=spectral_galerkin(Ax2,b2,s2,'adapt','solver',@pcg);
[X,errz]=spectral_galerkin(A2,b2,s2,'adapt','solver',@pcg,'lowmem',1);
[X,errz]=spectral_galerkin(Ax2,b2,s2,'adapt','solver',@pcg,'lowmem',1);
[X,errz]=spectral_galerkin(A2,b2,s2,'adapt','solver',@pcg,'lowmem',1,'ptol',1e-13);