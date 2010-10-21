function r = minimum_coefficient(X)
%MINIMUM_COEFFICIENT Type of error estimate for given 'X'
%
% r = minimum_coefficient(X)
% 
% Computes the average magnitude of the coefficients associated with the
% basis functions of the two highest degrees. 

% Copyright 2009-2010 David F. Gleich (dfgleic@sandia.gov) and Paul G. 
% Constantine (pconsta@sandia.gov)
%
% History
% -------
% :2010-06-14: Initial release


N=size(X.coefficients,1); dim=size(X.index_set,1);
maxorder=max(max(X.index_set));

if dim>1
    indz1=sum(X.index_set)==maxorder;
    indz2=sum(X.index_set)==(maxorder-1);
    r=(1/(2*N))*(sum(sum(X.coefficients(:,indz1).^2)) + ...
        sum(sum(X.coefficients(:,indz2).^2)));
else
    r=sum(sum(X.coefficients(:,end-1:end).^2))/(2*N);
end
r = sqrt(r);