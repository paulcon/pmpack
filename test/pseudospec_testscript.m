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
[X,errz]=pseudospectral(iAb,s,'adapt','ConvType','MinCoeff');
[X,errz]=pseudospectral(iAb2,s2,'adapt','ConvType','MinCoeff');
[X,errz]=pseudospectral(iAb,s,'adapt','ConvType','RelErr');
[X,errz]=pseudospectral(iAb2,s2,'adapt','ConvType','RelErr');
[X,errz]=pseudospectral(iAb,s,'adapt','RefSoln',Y);
[X,errz]=pseudospectral(iAb2,s2,'adapt','RefSoln',Y2);

[X,errz]=pseudospectral(iAb,s,'adapt','ConvType','resid','matfun',A,'vecfun',b);
[X,errz]=pseudospectral(iAb2,s2,'adapt','ConvType','resid','matfun',A2,'vecfun',b2);

[X,errz]=pseudospectral(iAb,s,'adapt','ConvType','resid','matvecfun',Ax,'vecfun',b);
[X,errz]=pseudospectral(iAb2,s2,'adapt','ConvType','resid','matvecfun',Ax2,'vecfun',b2);

matlabpool('open',2);
[X,errz]=pseudospectral(iAb,s,'adapt','ConvType','RelErr','parallel',1);
[X,errz]=pseudospectral(iAb2,s2,'adapt','ConvType','RelErr','parallel',1);
matlabpool('close');




