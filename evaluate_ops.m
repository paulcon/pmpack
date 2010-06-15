function P = evaluate_ops(s,n,points)
%EVALUATE_OPS Evaluate the n-vector of orthonormal polynomials at a point
%
% P = evaluate_ops(s,n,points) 
%
% Returns a vector of the first n orthonormal polynomials evaluated at 
% point or a set of points
%
% Inputs:
%   s     : parameter struct
%   n     : number of polynomials to evaluate.
%   point : the point at which to evaluate the vector of polynomials.
%
% Outputs:
%   P     : a vector of orthonormal polynomials evaluated at a point.
%
% Example:
%   % plot the first 5 Legendre polynomials
%   xx = [-1:0.1:1]; % generate a grid
%   P = evaluate_ops(parameter,5,xx); % evaluate on xx
%   plot(xx,P'); % plot
%
% See also EVALUATE_EXPANSION 


% Copyright 2009-2010 David F. Gleich (dfgleic@sandia.gov) and Paul G. 
% Constantine (pconsta@sandia.gov)
%
% History
% -------
% :2010-06-14: Initial release

assert(n>0, 'Please ensure n>0.');
assert(max(size(s))==1, 'Works for single parameters.'); 

P=zeros(n,length(points));
ab=s.recur(n);
P(1,:)=1/sqrt(ab(1,2));
if n==1; return; end
P(2,:)=(points(:)'-ab(1,1))./sqrt(ab(2,2));
if n==2; return; end
for i=3:n
    P(i,:)=((points(:)'-ab(i-1,1)).*P(i-1,:)-sqrt(ab(i-1,2))*P(i-2,:))./sqrt(ab(i,2));
end


