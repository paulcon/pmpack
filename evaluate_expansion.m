function [r,basis]=evaluate_expansion(X,p)
% EVALUATE_EXPANSION Evaluates a polynomial expansion at a point.
%
% r = evaluate_expansion(X,p);
% [r,basis] = evaluate_expansion(X,p);
%
% The input 'X' is a struct containing the variables for the polynomial
% approximation. The input 'p' is a point in the parameter space to
% evaluate the approximation.
%
% References:
%   Constantine, P.G., Gleich, D.F., Iaccarino, G. 'Spectral Methods for
%       Parameterized Matrix Equations'. arXiv:0904.2040v1, 2009.
%
% Copyright 2010 David F. Gleich (dfgleic@sandia.gov) and Paul G. 
% Constantine (pconsta@sandia.gov)

s=X.variables;

nd=size(X.index_set,1);          % the dimension of the point 
orders=max(X.index_set,[],2)'+1; % the highest expansion order for each dimension
% quick checks
assert(nd==length(p),'Please ensure nd=size(X.index_set,1) and length(point) are equal');
assert(nd==length(s),'Please ensure nd=length(X.variables) and length(point) are equal');

% Compute polynomials evaluated at each point.
polys = cellfun(@evaluate_ops,num2cell(s(:)),num2cell(orders(:)),num2cell(p(:)),...
    'UniformOutput',false);
% TODO Compare timing with the following loop code
% polys=cell(1,nd); for i=1:nd, polys{i}=evaluate_ops(x(i),orders(i),p(i)); end

nb=size(X.index_set,2);          % nb is the total number of basis polynomials
basis=ones(nb,1);        % basis is the polynomial basis evaluated at point
for i=1:nb
    for j=1:nd
        basis(i)=basis(i)*polys{j}(X.index_set(j,i)+1);
    end
end

r=X.coefficients*basis;             % compute the expansion via combinations of the basis
