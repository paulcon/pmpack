function Q=jacobi_eigenvecs(s,n)
% JACOBI_EIGENVECS Construct the eigenvectors of the Jacobi matrix.
%
% Q = jacobi_eigenvecs(s,n)
%
% Given a single parameter 's' and a scalar 'n', this function constructs
% the eigenvectors of the symmetric, tridiagonal matrix of recursion 
% coefficients for the orthogonal polynomial basis associated with the 
% measure on the space of 's'. 
%
% Example:
%   s = legendre_parameter();
%   Q = jacobi_eigenvecs(s,5);
%
% Copyright 2010 David F. Gleich (dfgleic@sandia.gov) and Paul G. 
% Constantine (pconsta@sandia.gov).


J=jacobi_matrix(s,n);
[Q,D]=eig(J);
[D,I]=sort(diag(D));
Q=Q(:,I);