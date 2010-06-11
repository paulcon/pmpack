function s = hermite_parameter()
% HERMITE_PARAMETER Construct a parameter with a Gaussian weight function.
%
% s = hermite_parameter() 
%
% Generates a parameter on the interval [-Inf,Inf] with a Gaussian weight 
% function.
%
% See also PARAMETER, JACOBI_PARAMETER, LEGENDRE_PARAMETER
%
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