function [X,errz] = spectral_galerkin(A,b,s,pOrder,varargin)
%SPECTRAL_GALERKIN Galerkin approximation of solution to A(s)x(s)=b(s)
%
% X = spectral_galerkin(A,b,s,pOrder);
% X = spectral_galerkin(A,b,s,pOrder,...);
% [X,err] = spectral_galerkin(A,b,s,pOrder,...);
%
% The function spectral_galerkin computes the Galerkin approximation to the
% solution x(s) of the parameterized matrix equation A(s)x(s)=b(s) using a
% basis of multivariate orthogonal polynomials. 
%
% Outputs:
%   X:          A struct containing the components of the Galerkin 
%               solution. See below for a more detailed description.
%
%   err:        An estimate of the error in the approximation. If pOrder 
%               is set to 'adapt', then this is a vector of the error 
%               estimates to examine convergence.
%
% Required inputs:
%   A:          A function handle that can take one of two forms. It
%               either returns the parameterized matrix evaluated at a
%               point, e.g. @(s) A(s). Or it returns the action of the
%               matrix evaluated at a point multiplied by a given vector,
%               e.g. @(s,v) A(s)*v.
%
%   b:          An anonymous function that returns the parameterized right
%               hand side evaluated at a point, e.g. @(s) b(s).
%
%   s:          A vector of parameter structs. The length of s is
%               considered to be the dimension d of the parameter space. 
%               See the function parameter.m.
%
%   pOrder:     The order of the polynomial approximation. A scalar input
%               creates a full polynomial basis set for the given 
%               dimension. A vector input creates a tensor product basis 
%               set with the order specified for each dimension by the 
%               components of pOrder. Set this to the string 'adapt' to 
%               increase the polynomial order of a full polynomial basis
%               until the chosen error estimate is below a given tolerance.
%               Finally, this can be an array of columns of nonnegative 
%               integers of length d representing the multi-indices that
%               define a multivariate basis polynonmial. See the examples
%               for more details.
%
% Optional inputs:
% To specify optional inputs, use the 'key/value' format. For example, to
% set the convergence tolerance 'pTol' to 1e-6, include 'pTol',1e-6 in 
% the argument list. See the examples below for more details.
%
%   qOrder:    Order of tensor product Gaussian quadrature integration 
%               used to construct the Galerkin system and the right hand 
%               side. Can be a scalar or a vector of positive integers of
%               length d. (Default 2*(pOrder+1))
%
%   Solver:     An anonymous function of the form @(A,b) solver(A,b) that
%               solves the constant matrix equation Ax=b. See the examples
%               for more details. (Default @mldivide)
%
%   LowMem:     A flag taking values 0 or 1 that determines how the solver
%               stores the matrices evaluated at the quadrature points. If
%               the input 'A' returns a matrix, then LowMem=1 tells the 
%               code to store each A(lambda) in a cell array, for the 
%               quadrature points lambda. Otherwise, the code explicitly
%               forms the Galerkin matrix and solves it with 'Solver'. If
%               the input 'A' is a matrix vector product interface, then
%               'LowMem' is ignored. (Default 0)
%
%   pTol:       A scalar representing the tolerance for the pOrder='adapt'
%               option. This is ignored if 'pOrder' is not set to 'adapt'.
%               (Default 1e-8)
%
%   ErrEst:     A string that determines the type of error estimate to use.
%               The options include: 'relerr' computes the difference 
%               between the approximation and a reference solution.
%               'mincoeff' computes the average of the magnitude of the
%               coefficients associated with the terms of the the two
%               highest degree. 'resid' uses the inputs 'A' and 'b' to
%               compute a residual error estimate. (Default 'relerr')
%
%   RefSoln:    A struct containing a reference solution to compare against
%               a computed approximation. If pOrder='adapt', then this is
%               set as the approximation with 'pOrder' one less than the
%               current approximation. (Default [])
%
%   ParRhs:     A flag taking values 0 or 1 that tells the code to
%               construct the Galerkin right hand side using the Parallel 
%               Computing Toolbox. (Default 0)
%
%   Verbose:    A flag taking values 0 or 1 that tells the code whether or
%               not to print detailed status information during the
%               computation of the approximation. (Default 0)
%
% The output struct 'X' contains the following fields.
%   X.coefficients: An array of size N by # of bases containing the
%               coefficients of the Galerkin approximation.
%
%   X.index_set: An array of size d by # of bases containing the
%               multi-indicies corresponding to each basis polynomial.
%
%   X.variables: The input vector of parameters 's' used to construct the
%               approximation.
%
%   X.fun:      If 'X' is a pseudospectral approximation, this is the
%               anonymous function used to compute the coefficients.
%
%   X.matfun:   If the input 'A' returns a matrix, this is a handle to that
%               function.
%
%   X.vecfun:   A handle to the input 'b'.
%
%   X.matvecfun: If the input 'A' returns a matrix-vector multiply, this is
%               a handle to that function.
%
% References:
%   Constantine, P.G., Gleich, D.F., Iaccarino, G. 'A factorization of 
%       the spectral Galerkin system for parameterized matrix equations: 
%       Derivation and applications'. arXiv: , 2010.
%
% Example:
%   A = @(t) [2 t; t 1];                    % 2x2 parameterized matrix
%   b = @(t) [2; 1];                        % constant right hand side
%   s = parameter();                        % parameter defined on [-1,1]
%   pOrder = 13;                            % degree 13 approximation
%   X = spectral_galerkin(A,b,s,pOrder);
%
% See also PSEUDOSPECTRAL


% Copyright 2009-2010 David F. Gleich (dfgleic@sandia.gov) and Paul G. 
% Constantine (pconsta@sandia.gov)
%
% History
% -------
% :2010-06-14: Initial release


if nargin<4, error('Not enough input arguments.'); end

% get constants
dim=length(s);
s_mid=midpoint(s); 
v=b(s_mid); 
N=length(v);

% find out what type of functions were passed.
e=[]; 
try A(s_mid,v); catch e, end
if isempty(e)
    matfun=[];
    matvecfun=A;
else
    if isequal(e.identifier,'MATLAB:TooManyInputs') 
        matfun=A;
        matvecfun=[]; 
    else
        error('Something went terribly wrong: %s',e.message);
    end
end
vecfun=b;

% set default values
qOrder=[];
solver=[];
lowmem=0;
pTol=0;
errest='relerr'; % types: relerr, mincoeff, resid
refsoln=[];
errz=[];
par_rhs=0;
verbose=0;
vprintf = @(varargin) fprintf('spectral_galerkin: %s\n',sprintf(varargin{:}));

for i=1:2:(nargin-4)
    switch lower(varargin{i})
        case 'qorder'
            qOrder=varargin{i+1};
            if ~isempty(qOrder)
                if isscalar(qOrder)
                    qOrder=qOrder*ones(1,dim);
                else
                    if length(qOrder)~=dim, error('Length of integration order must equal dimension'); end
                end
            end
        case 'solver'
            solver=varargin{i+1};        
        case 'lowmem'
            lowmem=varargin{i+1};
        case 'ptol'
            pTol=varargin{i+1};
        case 'errest'
            errest=lower(varargin{i+1});
        case 'refsoln'
            refsoln=varargin{i+1};
        case 'parrhs'
            par_rhs=varargin{i+1};
        case 'verbose'
            verbose=varargin{i+1};
        otherwise
            error('Unrecognized option: %s\n',varargin{i});
    end
end

if ~verbose, vprintf = @(varargin) []; end

% logic based on pOrder
if isnumeric(pOrder)
    if pTol~=0, warning('pmpack:optionIgnored','The specified polynomial tolerance will be ignored.'); end
    if isscalar(pOrder) % scalar order
        if isempty(qOrder), qOrder=2*(pOrder+1)*ones(dim,1); end
        basis=index_set('full',pOrder,dim);
    elseif min(size(pOrder))==1 % tensor order
        if max(size(pOrder))~=dim, error('Tensor order must equal dimension.'); end
        if isempty(qOrder), qOrder=2*(pOrder+1); end
        basis=index_set('tensor',pOrder);
    elseif size(pOrder,1)==dim % a given basis set
        basis=pOrder;
        pOrder=max(pOrder,[],2);
        if isempty(qOrder), qOrder=2*(pOrder+1); end
    else
        error('Unrecognized pOrder.');
    end
elseif isequal(pOrder,'adapt')
    if pTol==0, pTol=1e-8; end
else
    error('Unrecognized option for pOrder: %s\n',pOrder);
end



if isequal(pOrder,'adapt')
    vprintf('using adaptive computation ErrEst=%s',errest);
    
    if isequal(errest,'mincoeff') && ~isempty(refsoln)
        warning('pmpack:optionIgnored','Reference solution will be ignored.');
    end
    
    if isempty(refsoln)
        vprintf('computing reference solution');
        
        refsoln=spectral_galerkin(A,b,s,0,varargin{:});
    end
    
    err=inf; order=1;
    while err>pTol
        vprintf('adaptive solution order=%i, error=%g\n',order, err);
        X=spectral_galerkin(A,b,s,order,varargin{:});
        err=error_estimate(errest,X,refsoln);
        errz(order)=err; %#ok<AGROW>
        order=order+1;
        if isequal(errest,'relerr'), refsoln=X; end
    end
else
    if any(qOrder<(pOrder+1))
        qOrder=max(pOrder+1,qOrder);
        warning('pmpack:insufficientOrder',...
            ['Integration order insufficient, ' ...
             'qOrder must be larger than pOrder. Using pOrder instead.']);
    end
    
    vprintf('constructing quadrature rule npolys=%i max_qorder=%i',...
        prod(qOrder+1), max(qOrder));
    
    % Construct the gauss points and truncated Jacobi eigenvectors needed to
    % form the Galerkin system.
    Q=cell(dim,1);
    for i=1:dim
        Q{i}=jacobi_eigenvecs(s(i),qOrder(i));
    end 
    p=gaussian_quadrature(s,qOrder);

    vprintf('constructing generalized quadrature weights nbasis=%i',...
        size(basis,2));
    nbasis=size(basis,2);
    QQ=zeros(nbasis,size(p,1)); 
    for i=1:nbasis
        index=basis(:,i);
        qvec=1;
        for j=1:length(index)
            qvec=kron(qvec,Q{j}(index(j)+1,:));
        end
        QQ(i,:)=qvec;
    end
    
    % Construct the Galerkin right hand side.
    if par_rhs
        t0=tic;
        b0 = vecfun(p(1,:));
        dt=toc(t0);
        vprintf('constructing rhs with parfor npoints=%i time_for_one=%.1fs',...
            size(p,1),dt);
        BD = zeros(length(b0),size(p,1));
        BD(:,1) = b0*QQ(1,1);
        parfor pi=2:size(p,1)
            BD(:,pi) = vecfun(p(pi,:))*QQ(1,pi); %#ok<PFBNS>
        end
        Grhs=reshape(BD*QQ',nbasis*N,1);
    else
        B=cell2mat(cellfun(vecfun,num2cell(p,2),'UniformOutput',0)');
        D=spdiags(QQ(1,:)',0,size(QQ,2),size(QQ,2));
        Grhs=reshape((B*D)*QQ',nbasis*N,1);
    end
    
    if ~isempty(matfun)
        % matfun is specified (i.e. is "is not empty"), so 
        % we compute all the matrices once, and save them in a cell array.
        vprintf('precomputing matrices at sample points npoints=%i',...
            size(p,1));
        
        % on 2010-06-07, dgleich checked that this code was not generally
        % slower than using cellfun.
        Acell = cell(size(p,1),1);
        parfor pi=1:size(p,1)
            Acell{pi} = matfun(p(pi,:)); %#ok<PFBNS>
        end
        
        if lowmem
            Gfun=@(v) gmatvec_lowmem(v,Acell,QQ);
            
            if isempty(solver), solver=@(A,b)gmres(A,b); end
            vprintf('LowMem option matsize=%i solver=%s',...
                length(Grhs),func2str(solver));
            U=solver(Gfun,Grhs);
        else
            Alambda=blkdiag(Acell{:});
            Gmat=kron(QQ,speye(N))*Alambda*kron(QQ',speye(N));
            
            if isempty(solver), solver=@(A,b)mldivide(A,b); end
            vprintf('fullmat option matsize=%i nnz/row=%g solver=%s',...
                length(Grhs),nnz(Gmat)/size(Gmat,1),func2str(solver));
            U=solver(Gmat,Grhs);
        end
    elseif ~isempty(matvecfun)
        
        if lowmem 
            warning('pmpack:optionIgnored','Low memory option ignored.'); 
        end
        Gfun=@(v) gmatvec(v,matvecfun,QQ,p);
        
        if isempty(solver), solver=@(A,b)gmres(A,b); end
        vprintf('using matvecfun algorithms matsize=%i solver=%s',...
                length(Grhs),func2str(solver));
        U=solver(Gfun,Grhs);
    else
        error('Something went terribly wrong.');
    end
    vprintf('done');
    
    % Pack up the solution, dude.
    X.coefficients=reshape(U,N,nbasis);
    X.index_set=basis;
    X.variables=s; 
    X.fun=[];
    X.matfun=matfun;
    X.vecfun=vecfun;
    X.matvecfun=matvecfun;
    X=sort_bases(X);
end

end

function z=gmatvec_lowmem(v,Amats,Q)
% MATVEC Computes the product of the Galerkin matrix times a vector.
[m,n]=size(Q);
N=length(v)/m;

X=reshape(v,N,m);
X=X*Q;
refsoln=zeros(N,n);
parfor i=1:n
    refsoln(:,i)=Amats{i}*X(:,i);
end
Z=refsoln*Q';
z=Z(:);
end

function z=gmatvec(v,Ax,Q,p)
% MATVEC Computes the product of the Galerkin matrix times a vector.
[m,n]=size(Q);
N=length(v)/m;

X=reshape(v,N,m);
X=X*Q;
refsoln=zeros(N,n);
parfor i=1:n
    refsoln(:,i)=Ax(p(i,:),X(:,i)); %#ok<PFBNS>
end
Z=refsoln*Q';
z=Z(:);
end

