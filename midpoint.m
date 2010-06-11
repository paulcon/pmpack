function m = midpoint(s)
%MIDPOINT Computes the midpoint of the domain defined by 's'.

dim=length(s);
m=zeros(1,dim);
for i=1:dim
    if ~isfield(s,'r') || ~isfield(s,'l'), error('Cannot determine midpoint.'); end
    m(i)=0.5*(s(i).r-s(i).l)+s(i).l;
end
