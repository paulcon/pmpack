function Q=jacobi_eigenvecs(p,n)
% JACOBI_EIGENVECS Construct the eigenvector matrix of the Jacobi matrix 
% for a parameter.
%
% Given a parameter represented by a set of orthogonal polynomials, the
% Jacobi matrix J encodes the recursion coefficients to evaluate all the 
% orthogonal polynomials at a point. Since J is symmetric and tridiagonal,
% it admits a spectral decomposition J=QDQ' of real eigenvalues.
%
% Q = jacobi_eigenvecs(p,n) constructs the eigenvectors of the Jacobi 
% matrix of dimension n for parameter struct p.  
%
% Q = jacobi_eigenvecs(ab,n) constructs the eigenvectors of the Jacobi 
% matrix directly from the recursion coefficents.  This procedure requires 
% that size(ab,1)<=n.
%
% See also
%
% Example:
%   x = legendre_param();
%   Q = jacobi_eigenvecs(x,5);
%   
%
% Copyright, Stanford University, 2009
% Paul G. Constantine, David F. Gleich


J=jacobi_matrix(p,n);
[Q,D]=eig(J);
[D,I]=sort(diag(D));
Q=Q(:,I);