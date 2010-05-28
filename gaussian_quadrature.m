function [p,w]=gaussian_quadrature(s,n,gridflag)
% GAUSSRULE Compute a Gaussian quadrature rule
% 
% [p,w] = gaussian_quadrature(ab,n);
% [p,w] = gaussian_quadrature(jacobi_param(),n);
% [xs,ws] = gaussian_quadrature([jacobi_param(),jacobi_param()],n);
%
% See also
%
% Example:
%   
%
% Copyright, Stanford University, 2009
% Paul G. Constantine, David F. Gleich

if ~exist('gridflag','var') || isempty(gridflag), gridflag=1; end

if isstruct(s) 
    p = cell(size(s));
    w = cell(size(s));
    if isscalar(n), n=n*ones(size(s)); end
    for i=1:numel(s)
        [p{i},w{i}] = gaussian_quadrature(s(i).recur(n(i)),n(i));
    end

    if numel(p)==1,
        % dereference the cell if there is only one element
        p=p{1};
        w=w{1};
    else
        if gridflag
            % build a tensor grid of points
            pgrid=1; wgrid=1;
            for i=1:numel(p)
                pgrid=[kron(pgrid,ones(length(p{i}),1)) kron(ones(size(pgrid,1),1),p{i})];
                wgrid=kron(wgrid,w{i});
            end
            p=pgrid(:,2:end);
            w=wgrid;
        end
    end
else
    % s is a set of recursion coefficients
    if n>size(s,1)
        error('gaussrule:insufficientCoefficients',...
            ['The gaussrule needs more coefficients than points,' ...
                'but numcoeff=%i and numpoints (n)=%i'], size(s,1), n);
    end
        
    ab=s;
    J=jacobi_matrix(s,n); 
    [V,D]=eig(J);
    [p,I]=sort(diag(D));
    w=V(1,I)'.^2;
    w=ab(1,2)*w;
end
        
