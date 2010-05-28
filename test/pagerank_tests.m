%% PageRank test cases borrowed from rapr codes
msgid = 'pmpack:pagerankTests';

%% First check a known solution to make sure it's correct
P = [0 1/2 1/2;
     0 0 1;
     0 0 1];
xfun = @(a) (eye(3)-a*P')\((1-a)/3*[1;1;1]);
Afun = @(a) (eye(3)-a*P');
bfun = @(a) ((1-a)/3*[1;1;1]);
xtrue = @(a) [(1-a)/3; (2-a-a^2)/6; (2+3*a+a^2)/6];
for a = linspace(0,1-10*eps(1),50)
    if norm(xfun(a)-xtrue(a),'inf')>2*eps(1)
        error(msgid,'incorrect coded solution');
    end
end

%% Check against pseudo-spectral solutions
X = pseudospectral(xfun,jacobi_parameter(0,1,0,0),2);
for a = linspace(0,1-10*eps(1),50)
    if norm(evaluate_expansion(X,a)-xtrue(a),'inf')>10*eps(1)
        error(msgid,'incorrect pseudo-spectral solution');
    end
end

%% Check against a spectral galerkin solution
X = spectral_galerkin(Afun,bfun,jacobi_parameter(0,1,0,0),2);
for a = linspace(0,1-10*eps(1),50)
    if norm(evaluate_expansion(X,a)-xtrue(a),'inf')>10*eps(1)
        error(msgid,'incorrect spectral-galerkin solution');
    end
end

%% Compare spectral_galerkin to pseudo-spectral
% For PageRank, spectral_galerkin and pseudospectral produce the same
% expansion.  Let's verify our codes have this property
load('../demo/wb-cs.stanford.mat');
P = Pcc; % only use the largest strong component -- this fixes a few technical details
Afun = @(a) (speye(size(P)) - a*P');
bfun = @(a) (1-a)./size(P,1)*ones(size(P,1),1);
xfun = @(a) Afun(a)\bfun(a);

% this solution takes about 11 polynomials to resolve
for i=1:3:15
    X = pseudospectral(xfun,jacobi_parameter(0,1,2,3),i);
    Y = spectral_galerkin(Afun,bfun,jacobi_parameter(0,1,2,3),i);
    if norm(X.coefficients - Y.coefficients,'fro')>100*eps(1)
        error(msgid, ['spectral-galerkin and pseudo-spectral solutions do not match. ' ...
            'difference = %e'], norm(X.coefficients - Y.coefficients,'fro'));
    end
end
