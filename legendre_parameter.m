function s = legendre_parameter(l,r)
% LEGENDRE_PARAMETER Construct a parameter with a uniform weight function.
%
% s = legendre_parameter() 
% s = legendre_parameter(l,r) 
%
% Generates a parameter on the interval [l,r] with a uniform weight 
% function.
%
% See also PARAMETER, JACOBI_PARAMETER, HERMITE_PARAMETER
%
% Copyright 2010 David F. Gleich (dfgleic@sandia.gov) and Paul G. 
% Constantine (pconsta@sandia.gov).

s = jacobi_parameter(l,r,0,0);

end
