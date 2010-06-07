function [X,errz] = spectral_galerkin(A,b,s,p_order,varargin)

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
qorder=[];
solver=[];
lowmem=0;
ptol=0;
convtype='relerr'; % types: relerr, mincoeff, resid
refsoln=[];
errz=[];
parallel_rhs=0;
verbose=0;
vprintf = @(varargin) fprintf('spectral_galerkin: %s\n',sprintf(varargin{:}));

for i=1:2:(nargin-4)
    switch lower(varargin{i})
        case 'qorder'
            qorder=varargin{i+1};
            if ~isempty(qorder)
                if isscalar(qorder)
                    qorder=qorder*ones(1,dim);
                else
                    if length(qorder)~=dim, error('Length of integration order must equal dimension'); end
                end
            end
        case 'solver'
            solver=varargin{i+1};        
        case 'lowmem'
            lowmem=varargin{i+1};
        case 'ptol'
            ptol=varargin{i+1};
        case 'convtype'
            convtype=lower(varargin{i+1});
        case 'refsoln'
            refsoln=varargin{i+1};
        case 'parallel_rhs'
            parallel_rhs=varargin{i+1};
        case 'verbose'
            verbose=varargin{i+1};
        otherwise
            error('Unrecognized option: %s\n',varargin{i});
    end
end

if ~verbose, vprintf = @(varargin) []; end

% logic based on p_order
if isnumeric(p_order)
    if ptol~=0, warning('The specified polynomial tolerance will be ignored.'); end
    if isscalar(p_order) % scalar order
        if isempty(qorder), qorder=2*(p_order+1)*ones(dim,1); end
        basis=index_set('full',p_order,dim);
    elseif min(size(p_order))==1 % tensor order
        if max(size(p_order))~=dim, error('Tensor order must equal dimension.'); end
        if isempty(qorder), qorder=2*(p_order+1); end
        basis=index_set('tensor',p_order);
    elseif size(p_order,1)==dim % a given basis set
        if isempty(qorder), qorder=2*(max(p_order,[],2)+1); end
        basis=p_order;
    else
        error('Unrecognized p_order.');
    end
elseif isequal(p_order,'adapt')
    if ptol==0, ptol=1e-8; end
else
    error('Unrecognized option for p_order: %s\n',p_order);
end



if isequal(p_order,'adapt')
    vprintf('using adaptive computation convtype=%s',convtype);
    
    if isequal(convtype,'mincoeff') && ~isempty(refsoln)
        warning('Reference solution will be ignored.');
    end
    
    if isempty(refsoln)
        vprintf('computing reference solution');
        
        refsoln=spectral_galerkin(A,b,s,0,...
            'qorder',qorder,'solver',solver,'lowmem',lowmem);
    end
    
    err=inf; order=1;
    while err>ptol
        vprintf('adaptive solution order=%i, error=%g\n',order, err);
        X=spectral_galerkin(A,b,s,order,varargin{:});
        err=error_estimate(convtype,X,refsoln);
        errz(order)=err;
        order=order+1;
        if isequal(convtype,'relerr'), refsoln=X; end
    end
else
    if any(qorder<(p_order+1))
        qorder=max(p_order+1,qorder);
        warning('pmpack:insufficientOrder',...
            ['Integration order insufficient, ' ...
             'qorder must be larger than p_order. Using p_order instead.']);
    end
    
    vprintf('constructing quadrature rule npolys=%i max_qorder=%i',...
        prod(qorder+1), max(qorder));
    
    % Construct the gauss points and truncated Jacobi eigenvectors needed to
    % form the Galerkin system.
    Q=cell(dim,1);
    for i=1:dim
        Q{i}=jacobi_eigenvecs(s(i),qorder(i));
    end 
    p=gaussian_quadrature(s,qorder);

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
    if parallel_rhs
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
            vprintf('lowmem option matsize=%i solver=%s',...
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
    
    % Pack up the solution.
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
for i=1:n
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
for i=1:n
    refsoln(:,i)=Ax(p(i,:),X(:,i));
end
Z=refsoln*Q';
z=Z(:);
end

