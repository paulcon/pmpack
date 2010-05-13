function r = minimum_coefficient(X)
N=size(X.coefficients,1); dim=size(X.index_set,1);
if dim>1
    indz1=find(sum(X.index_set)==max(sum(X.index_set)));
    indz2=find(sum(X.index_set)==max(sum(X.index_set))-1);
    r=(1/(2*N))*(sum(sum(X.coefficients(:,indz1).^2)) + ...
        sum(sum(X.coefficients(:,indz2).^2)));
else
    r=sum(sum(X.coefficients(:,end-1:end).^2))/(2*N);
end