% utils test script

s=parameter();
ab=s.recur(5);
s=parameter('legendre');
ab=s.recur(5);
s=parameter('jacobi');
ab=s.recur(5);
s=parameter('jacobi',0,1,2,2);
ab=s.recur(5);
s=parameter('chebyshev');
ab=s.recur(5);
s=parameter('hermite');
ab=s.recur(5);

J=jacobi_matrix(s,5);
Q=jacobi_eigenvecs(s,5);

s=[parameter(); parameter(); parameter()];

[p,w]=gaussian_quadrature(s,5);
[p,w]=gaussian_quadrature(s(1),5);
[p,w]=gaussian_quadrature(s,[5 4 3]);
[p,w]=gaussian_quadrature(s,[5 4 3],0);
[p,w]=gaussian_quadrature(s,[5 4 3],1);


I=index_set('tensor',4);
I=index_set('tensor',[4 3 4]);
I=index_set('full',4,3);
I=index_set('total order',4,3);
I=index_set('complete',4,3);
c=@(a)sum(a)<=4; I=index_set('constrained',4,3,c);
c=@(a)norm(a,inf)<=4; I=index_set('constrained',4,3,c);
w=[1 2]; c=@(a) a(:)'*w(:)<=6; I=index_set('constrained',6,2,c);

