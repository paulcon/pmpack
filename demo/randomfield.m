function [F,KL] = randomfield(cv,mesh,varargin)
%RANDOMFIELD Generates realizations of a Gaussian random field.
%
% F = randomfield(cv,mesh)
% F = randomfield(cv,mesh,...)
% [F,KL] = randomfield(cv,mesh,...)
%   
% Outputs:
%   F:      A matrix of random field realizations. Each column is a
%           realization of the random field with each element corresponding
%           to a point supplied in the required input 'mesh'.
%
%   KL:     A struct containing the components of a Karhunen-Loeve (KL)
%           expansion of the random field. See below for a more detailed
%           description.
%
% Required inputs:
%   cv:         Either an anonymous function representing a two point
%               covariance function of the form cv(x1,x2) or a covariance
%               matrix associated with the mesh. See 'cvfun' and 'cvmat'
%               below for more details.
%
%   mesh:       A matrix of size nx by d, where nx is the number of points 
%               in the mesh and d is the dimension.
%
% Optional inputs:
% To specify optional inputs, use the 'key/value' format. For example, to
% set the number of realizations to 7, include 'nsamples',7 in the argument
% list. See the examples below for more details.
%
%   nsamples:   Number of realizations of the random field. (Default 1)
%
%   data:       A struct with prespecified data values that yields a
%               conditional random field. See below for the details of the
%               data struct. (Default [])
%
%   filter:     Set filter to a number between 0 and 1 to capture a
%               percentage of the energy of the field as determined by the
%               singular values of the covariance matrix. (Default 1 means 
%               no filter)
%
%   trunc:      Explicitly sets a truncation level for the KL expansion.
%               (Default 0 means no truncation)
%
%   spthresh:   Threshold for setting an element of the covariance matrix
%               to zero in the sparse matrix. (Default 0 means include all
%               elements).
%
%   mean:       A user supplied mean of the random field. Must match the
%               values of the given mesh. (Default zeros)
%
%   cvfun:      A function handle to a two point covariance function. The
%               function handle must be of the form cv=@(x1,x2) 
%               my_cov(x1,x2), where x1 and x2 are row vectors giving the 
%               coordinates of two points in the given mesh. (Default [])
%
%   cvmat:      A struct of user supplied covariance matrices. This is 
%               useful for large meshes if one wishes to form the matrix 
%               once and store it for quick retrieval. See below for 
%               details on the fields in the struct. (Default [])
%
%   snaps:      A data matrix of realizations (snapshots) from a random
%               field on the mesh, where each column is a realization.
%
% The output struct 'KL' with components of the Karhunen-Loeve representation
% of the random field has the following fields.
%   KL.mean:    The mean of the random field. If mean was supplied by the
%               user, then this is the same vector. If the field was
%               conditioned on data, then this is equivalent to Kriging
%               interpolant.
%
%   KL.bases:   The eigenvectors covariance matrix.
%
%   KL.sv:      The square root of the eigenvalues of the covariance
%               matrix. Plot these on a semilog scale to examine their
%               decay.
%
% The input struct 'data' must contain the following fields.
%   data.x:     A matrix of size ndata by d, where ndata is the number of
%               data points and d is the dimension of the mesh, containing
%               the points in the mesh where the field is known.
%
%   data.fx:    A vector of length ndata containing the known values of 
%               the random field.
%
% The input struct 'cvmat' may contain any of the following fields.
%   cvmat.C:    The covariance matrix between the unknown elements of the 
%               field.
%
%   cvmat.A:    The covariance matrix between data points.
%
%   cvmat.B:    The covariance matrix between data points and unknowns. The
%               code expects this to be structured so that rows correspond
%               to mesh points and columns correspond to data points.
%
% Example:
%   cv = @(x1,x2) gp_exp_cov(x1,x2,1,1);    % exponential covariance
%   mesh = linspace(100,-1,1);              % generate a mesh
%   data.x = [-1; 1]; data.fx = [0; -1];    % specify boundaries
%
%   [F,KL] = randomfield(cv, mesh, ...
%               'nsamples', 10, ...
%               'data', data, ...
%               'filter', 0.95);
%
%   % to generate 100 more samples using the KL
%   trunc = length(KL.sv);                  % get the truncation level
%   W = randn(trunc,100); 
%   F2 = KL.mean*ones(1,100) + KL.bases*KL.sv*W;
%   
% References:
%   M. Davis, "Production of Conditional Simulations via the LU Triangular
%       Decomposition of the Covariance Matrix". Mathematical Geology,
%       1987.
%
% Copyright 2010 Qiqi Wang (qiqi@mit.edu) and Paul G. Constantine 
% (pconsta@sandia.gov).

if nargin<2, error('Not enough inputs.'); end

nx=size(mesh,1);

% set default values.
nsamples=1;
data=[];
filter=1;
trunc=0;
spthresh=0;
mu=zeros(nx,1);
cvfun=[];
cvmat=struct([]);
X=[];

% determine the type of covariance information.
if isnumeric(cv)
    cvmat(1).C=cv;
elseif isstruct(cv)
    cvmat=cv;
else
    cvfun=cv;
end

% set input values and some error checking
for i=1:2:(nargin-2)
    switch lower(varargin{i})
        case 'nsamples'
            nsamples=varargin{i+1};
            if nsamples<1, error('nsamples must be greater than 1.'); end
        case 'data'
            data=varargin{i+1};
            if ~isfield(data,'x'), error('data struct must have field x.'); end
            if ~isfield(data,'fx'), error('data struct must have field fx.'); end
            if size(mesh,2) ~= size(data.x,2)
                error('Data and mesh are different dimensions.');
            end
        case 'filter'
            filter=varargin{i+1};
            if filter<=0 || filter>1, error('filter must be strictly between 0 and 1.'); end
        case 'trunc'
            trunc=varargin{i+1};
            if trunc<0 || trunc>nx, error('trunc must be positive and smaller than mesh size.'); end
        case 'spthresh'
            spthresh=varargin{i+1};
            if spthresh<0, error('spthresh must be positive.'); end
        case 'mean'
            mu=varargin{i+1};
            if size(mu,1)~=nx, error('mean must match mesh.'); end
        case 'cvfun'
            cvfun=varargin{i+1};
        case 'cvmat'
            cvmat=varargin{i+1};
            if isfield(cvmat,'C')
                if size(cvmat.C,1)~=nx, error('Covariance matrix must match mesh.'); end
            end
        case 'snaps'
            X=varargin{i+1};
            if size(X,1)~=nx, error('snaps must match mesh.'); end
        otherwise
            error('Unrecognized option: %s\n',varargin{i});
    end
end

% form an empirical covariance matrix from data and set a truncation level
if ~isempty(X)
    m=size(X,2);
    cvmat(1).C=(1/(m-1))*(X*X');
    if trunc==0
        trunc=m-1;
    else
        trunc=min(trunc,m-1);
    end
end

if isempty(cvfun) && ~isfield(cvmat,'C')
    error('No covariance information is provided for the given mesh.');
end

% infer an exponential covariance model from the given covariance matrix
if isempty(cvfun) && isfield(cvmat,'C')
    [c,sigma]=infer_correlation(mesh,cvmat.C);
    cvfun=@(x1,x2) gp_exp_cov(x1,x2,c,sigma);
end

% construct covariance matrix between unknowns
if isfield(cvmat,'C');
    C=sparse(cvmat.C);
else
    C=covariance_matrix(cvfun,mesh,[],spthresh);
end

if ~isempty(data)
    if isfield(cvmat,'A');
        A=sparse(cvmat.A);
        if size(A,1)~=size(data.x,1)
            error('Data-data covariance matrix does not match data mesh.');
        end
    else
        A=covariance_matrix(cvfun,data.x,[],spthresh);
    end
    
    if isfield(cvmat,'B');
        B=sparse(cvmat.B);
        if size(B,1)~=size(mesh,1)
            error('Data-unknowns covariance matrix does not match mesh.');
        end
        if size(B,2)~=size(data.x,1)
            error('Data-unknowns covariance matrix does not match data mesh.');
        end
    else
        B=covariance_matrix(cvfun,mesh,data.x,spthresh);
    end

    % build the tranformation
    L=chol(A,'lower');
    Bpt=L\(B'); Bp=Bpt';
    C=C-Bp*Bpt;
    
    % compute the mean
    mu=Bp*(L\data.fx);
end

% compute the transform
if trunc
    [U,S,V]=svds(C,trunc);
elseif filter<1
    s=svd(full(C)); ss=s/sum(s);
    trunc=sum(cumsum(ss)<filter);
    [U,S,V]=svds(C,trunc);
else
    [U,S,V]=svd(full(C));
    trunc=nx;
end
T=U*diag(sqrt(diag(S)));

% generate the samples
W=randn(trunc,nsamples);
F=T*W+mu*ones(1,nsamples);

% construct the KL representation
KL.mean=mu;
KL.bases=U; 
KL.sv=sqrt(diag(S));




