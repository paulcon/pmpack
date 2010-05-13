function s = jacobi_parameter(l,r,a,b)
% JACOBI_PARAMETER Construct a parameter with a Jacobi weight function
%
% x = jacobi_param() generates a new Jacobi parameter over [-1,1] with 
% a uniform weight.  This parameter is the canonical integration parameter.
%
% x = jacobi_param(l,r,a,b) generate a new Jacobi parameter over [l,r]
% with parameters alpha=a and beta=b.
%
% See also LEGENDRE_PARAM
%
% Example:
%   % Integrate exp(-x^2) over -1,1
%   x = jacobi_param();
%   xw = gaussrule(x,5); % generate a 5 point Gaussian quadrature 
%   I = sum(feval(@(x) exp(-x.^2), xw(:,1)).*xw(:,2))
%
% Copyright, Stanford University, 2009
% Paul G. Constantine, David F. Gleich

if ~exist('l','var') || isempty(l), l= -1; end
if ~exist('r','var') || isempty(r), r=  1; end
if ~exist('a','var') || isempty(a), a=  0; end
if ~exist('b','var') || isempty(b), b=  0; end

s.name = sprintf('jacobi(%g,%g) with support [%g,%g]', a, b, l, r);
s.recur = @(n) jacobi_recur(n,l,r,a,b);
s.l=l; s.r=r;
end

function ab=jacobi_recur(n,l,r,a,b)
% Compute the recurrence coefficients for the Jacobi polynomials
a0 = (b-a)/(a+b+2);
ab = zeros(n,2);
b2a2 = b^2 - a^2;
s = (r-l)/2; o = l + (r-l)/2;
if n>0
    ab(1,1) = s*a0+o;
    ab(1,2) = 1;
end
for k=2:n
    ab(k,1) = s*b2a2/((2*(k-1)+a+b)*(2*k + a+b))+o;
    ab(k,2) = ((r-l)^2*(k-1)*(k-1+a)*(k-1+b)*(k-1+a+b)) / ...
                    ((2*(k-1)+a+b)^2*(2*(k-1)+a+b+1)*(2*(k-1)+a+b-1));
end

end
