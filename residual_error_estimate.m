function r=residual_error_estimate(X,Av,bfun)
%RESIDUAL_ERROR_ESTIMATE Computes residual error estimate
%
% r = residual_error_estimate(X)
% r = residual_error_estimate(X,Av,b)
%
% Computes the residual error estimate of the approximation 'X' using the
% operators A(s) and b(s). 
%
% See also ERROR_ESTIMATE MINIMUM_COEFFICIENT

% Copyright 2009-2010 David F. Gleich (dfgleic@sandia.gov) and Paul G. 
% Constantine (pconsta@sandia.gov)
%
% History
% -------
% :2010-06-14: Initial release

if nargin==1
    if ~isempty(X.matvecfun)
        Ax=X.matvecfun;
    elseif ~isempty(X.matfun)
        Ax=@(p,x) X.matfun(p)*x;
    else
        error('Unable to compute residual error estimate: Missing matvec or matvecfun');
    end
else
    Ax = Av;
end

if nargin==1
    if ~isempty(X.vecfun)
        b=X.vecfun;
    else
        error('Unable to compute residual error estimate: Missing vecfun.');
    end
else
    b = bfun;
end

% determine the order of Gauss quadrature necessary to evaluate the
% integral of the residual.
expansion_order=max(X.index_set,[],2);
gauss_order=4*(expansion_order+1);
[p,w]=gaussian_quadrature(X.variables,gauss_order);

r=0;
parfor i=1:size(p,1)
    u=evaluate_expansion(X,p(i,:));
    r=r+norm(Ax(p(i,:),u)-b(p(i,:)),2)^2*w(i); %#ok<PFBNS>
end

r=sqrt(r)/size(X.coefficients,1);
