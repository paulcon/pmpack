function s = parameter(type,l,r,a,b)
% PARAMETER Construct a parameter
%
% Copyright, Stanford University, 2010
% Paul G. Constantine, David F. Gleich

if ~exist('l','var') || isempty(l), l= -1; end
if ~exist('r','var') || isempty(r), r=  1; end
if ~exist('a','var') || isempty(a), a=  0; end
if ~exist('b','var') || isempty(b), b=  0; end

if nargin==0, type='legendre'; end

switch lower(type)
    case 'jacobi'
        s = jacobi_parameter(l,r,a,b);
    case 'legendre'
        s = jacobi_parameter(l,r,0,0);
    case 'chebyshev'
        s = jacobi_parameter(l,r,-1/2,-1/2);
    case 'hermite'
        s = hermite_parameter();
    otherwise
        error('Unrecognized parameter type: %s',type);
end

end
