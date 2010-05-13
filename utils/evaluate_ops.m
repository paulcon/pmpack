function P = evaluate_ops(s,n,point)
% EVALUATE_OPS Evaluate the n-vector of orthonormal polynomials at a point.
%
% P = evaluate_ops(x,n,point) returns a vector of the first n orthonormal
% polynomials evaluated at point.
%
% Inputs:
%   s     : parameter struct
%   n     : number of polynomials to evaluate.
%   point : the point at which to evaluate the vector of polynomials.
%
% Outputs:
%   P     : a vector of orthonormal polynomials evaluated at a point.
% 
% See also
%
% Example:
%   
%
% Copyright, Stanford University, 2009
% Paul G. Constantine, David F. Gleich

assert(n>0, 'Please ensure n>0.');
assert(max(size(s))==1, 'Works for single parameters.'); 

P=zeros(n,1);
ab=s.recur(n);
P(1)=1/sqrt(ab(1,2));
if n==1; return; end
P(2)=(point-ab(1,1))/sqrt(ab(2,2));
if n==2; return; end
for i=3:n
    P(i)=((point-ab(i-1,1))*P(i-1)-sqrt(ab(i-1,2))*P(i-2))/sqrt(ab(i,2));
end


