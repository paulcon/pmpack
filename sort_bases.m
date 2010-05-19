function Y=sort_bases(X)
% SORT_BASES Sorts the multivariate bases in order of total degree of basis
% elements.
%
% Y = sort_bases(X) sorts the basis elements of the expansion X by total 
% order and order the coefficients to match.
%
% See also 
%
% Example:
%   
%
% Copyright, Stanford University, 2009
% Paul G. Constantine, David F. Gleich

s=sum(X.index_set,1);
[t,ind]=sort(s,'ascend');
X.coefficients=X.coefficients(:,ind); 
X.index_set=X.index_set(:,ind);
Y=X;