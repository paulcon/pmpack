function s = hermite_parameter()
%HERMITE_PARAMETER Construct a parameter with a Gaussian weight function
%
% s = hermite_parameter() 
%
% Generates a parameter on the interval [-Inf,Inf] with a Gaussian weight 
% function.  This weight is slightly different from the standard Gaussian.
% It has variance 1/2 instead of 1.  
%
% Example: 
%   % see that we have variance 1/2
%   [p,w] = gaussian_quadrature(s,2);
%   std2=0; for i=1:length(p), std2=std2+p(i).^2*w(i); end, std2
%
% See also PARAMETER, JACOBI_PARAMETER, LEGENDRE_PARAMETER

% Copyright 2010 David F. Gleich (dfgleic@sandia.gov) and Paul G. 
% Constantine (pconsta@sandia.gov).

s.name = sprintf('hermite');
s.recur = @(n) hermite_recur(n);

end

function ab=hermite_recur(n)
% Compute the recurrence coefficients for the Hermite polynomials using
% Gautschi's OPQ software.
ab=r_hermite(n);
ab(1,2)=1;
end
