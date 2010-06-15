function m = midpoint(s)
%MIDPOINT Computes the midpoint of the domain defined by 's'
%
% s0 = midpoint(s)
%
% Sets s0 = [midpoint(s1) midpoint(s2) ... midpoint(sd)], that is,
% the midpoint of each dimension of your parameter space.

% Copyright 2009-2010 David F. Gleich (dfgleic@sandia.gov) and Paul G. 
% Constantine (pconsta@sandia.gov)
%
% History
% -------
% :2010-06-14: Initial release


dim=length(s);
m=zeros(1,dim);
for i=1:dim
    if ~isfield(s,'r') || ~isfield(s,'l'), error('Cannot determine midpoint.'); end
    m(i)=0.5*(s(i).r-s(i).l)+s(i).l;
end
