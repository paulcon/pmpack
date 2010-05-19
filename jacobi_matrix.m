function J=jacobi_matrix(p,n)
% JACOBI_MATRIX Construct a Jacobi matrix for a parameter.
%
% Given a parameter represented by a set of orthogonal polynomials, the
% Jacobi matrix J encodes the recursion coefficients to evaluate all the 
% orthogonal polynomials at a point.
%
% J = jacobi_matrix(p,n) constructs the Jacobi matrix of dimension n for
% parameter p.  
%
% J = jacobi_matrix(ab,n) constructs the Jacobi matrix directly from the
% recursion coefficents.  This procedure requires that size(ab,1)<=n.
%
% See also
%
% Example:
%   x = legendre_param();
%   J = jacobi_matrix(x,5);
%
% Copyright, Stanford University, 2009
% Paul G. Constantine, David F. Gleich

if ~isstruct(p) && (~exist('n','var') || isempty(n)), n=size(p,1); end
if isstruct(p)
    ab = p.recur(n);
elseif size(p,2)==2 && isfloat(p)
    ab = p;
    if n>size(ab,1)
        error('jacobi_matrix:insufficientCoefficients', ...
            ['Please increase the number of recursion coefficients' ...
                ' (currently %i) to construct the matrix of size %i'], ...
                size(ab,1), n);
    end
else 
    error('jacobi_matrix:invalidArgument', ...
        ['The parameter or set of recursion coefficients must be a struct',...
            ' or a n-by-2 matrix.']);
end
J = zeros(n,n);
% special case for n=1
if n==1, J = ab(1,1); return; end

J(1,1)=ab(1,1);
J(1,2)=sqrt(ab(2,2));
for i=2:n-1
    J(i,i)=ab(i,1);
    J(i,i-1)=sqrt(ab(i,2));
    J(i,i+1)=sqrt(ab(i+1,2));
end
J(n,n)=ab(n,1);
J(n,n-1)=sqrt(ab(n,2));
