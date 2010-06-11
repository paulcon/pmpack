% test script for pseudospectral method


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


[X,errz]=pseudospectral(iAb,s,5);
[X,errz]=pseudospectral(iAb2,s2,[5 5]);
[X,errz]=pseudospectral(iAb2,s2,[5 7]);
[X,errz]=pseudospectral(iAb,s,'adapt');
[X,errz]=pseudospectral(iAb,s,'adapt','ErrEst','MinCoeff');
[X,errz]=pseudospectral(iAb2,s2,'adapt','ErrEst','MinCoeff');
[X,errz]=pseudospectral(iAb,s,'adapt','ErrEst','RelErr');
[X,errz]=pseudospectral(iAb2,s2,'adapt','ErrEst','RelErr');
[X,errz]=pseudospectral(iAb,s,'adapt','RefSoln',Y);
[X,errz]=pseudospectral(iAb2,s2,'adapt','RefSoln',Y2);

[X,errz]=pseudospectral(iAb,s,'adapt','ErrEst','resid','matfun',A,'vecfun',b);
[X,errz]=pseudospectral(iAb2,s2,'adapt','ErrEst','resid','matfun',A2,'vecfun',b2);

[X,errz]=pseudospectral(iAb,s,'adapt','ErrEst','resid','matvecfun',Ax,'vecfun',b);
[X,errz]=pseudospectral(iAb2,s2,'adapt','ErrEst','resid','matvecfun',Ax2,'vecfun',b2);




