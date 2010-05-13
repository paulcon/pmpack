function [r,basis]=evaluate_expansion(X,p)
% EVALUATE_EXPANSION Evaluates a polynomial expansion at a point.
%
% See also
%
% Example:
%   
%
% Copyright, Stanford University, 2009
% Paul G. Constantine, David F. Gleich

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
