function [r,basis]=evaluate_expansion(X,p)
%EVALUATE_EXPANSION Evaluates a polynomial expansion at a point
%
% r = evaluate_expansion(X,p);
% [r,basis] = evaluate_expansion(X,p);
%
% The input 'X' is a struct containing the variables for the polynomial
% approximation. The input 'p' is a point in the parameter space to
% evaluate the approximation.
%
% Example:
%   % display a response surface
%   P = pmpack_problem('twobytwo');
%   X = pseudospectral(P.solve,P.s,5);
%   pts = [-1:0.01:1]; fvals_est=zeros(2,length(pts)); fvals_true=fvals_est;
%   for i=1:length(pts)
%     fvals_true(:,i) = P.solve(pts(i));
%     fvals_est(:,i) = evaluate_expansion(X,pts(i));
%   end
%   plot(pts,fvals_true,'k-',pts,fvals_est,'r-');
%
% See also ERROR_ESTIMATE EVALUATE_OPS

% Copyright 2009-2010 David F. Gleich (dfgleic@sandia.gov) and Paul G. 
% Constantine (pconsta@sandia.gov)
%
% History
% -------
% :2010-06-14: Initial release


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
