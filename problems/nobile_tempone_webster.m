function P = nobile_tempone_webster(varargin)
% NOBILE_TEMPONE_WEBSTER Construct an example from Nobile et al.
%
% This function requies the PDE toolbox.  It also produces a warning when
% we generate the mesh, you can safely ignore this warning.
%
% P = nobile_tempone_webster(); constructs a problem with 5 parameters
% P = nobile_tempone_webster(d); constructs a problem with d parameters
% P = nobile_tempone_webster('dim',d); constructs a problem with d
% parameters
% P = nobile_tempone_websiter('dim',d,'
%
% Return:
% The output P is a struct with the following fields:
% P.A - a parameterized matrix interface
% P.Av - a parameterized matrix-vector interface
% P.b - a parameterized right hand side
% P.N - the size of the matrix
% P.d - the number of parameters (1) for this problem
% P.solve - a parameterized solution
% P.s - the vector of parameters for this problem
% P.mesh - details on the mesh (points, edges, triangles);
%
% Example:
%   P = nobile_tempone_webster();
%   [X,r] = spectral_galerkin(P.Av,P.b,P.s,'adapt');
%   scatter(P.mesh.p(:,1),P.mesh.p(:,2),1,X.coefficient(:,1));
%

% Copyright, Stanford University, 2009-2010
% Paul G. Constantine, David F. Gleich

d = [];
if nargin==1,
    d = varargin{1};
else
    args = struct(varargin{:});
    if isfield(args,'dim'), d=args.dim; dim.size = []; end
end
if isempty(d), d=5; end

% Geometry and boundary conditions
g=[2 2 2 2;0 1 1 0;1 1 0 0;0 0 1 1;0 1 1 0;1 1 1 1;0 0 0 0];
b=[1 1 1 1;1 1 1 1;1 1 1 1;1 1 1 1;1 1 1 1;1 1 1 1;48 48 48 48;
    48 48 48 48;49 49 49 49;48 48 48 48];
[p,e,t]=initmesh(g,'Hmax',0.03);

% set these as nested functions so they share all the loaded data
    function [K,F]=ntw_assemble(s) 
        if length(s)~=d, error('pmpack:wrongParameterSize',...
                'this ntw problem takes %i parameters, not %i', d, length(s)); end

        % Set up coefficients.
        c=@(p,t,u,time) ntw_coeff(p,t,s);

        % Solve the system.
        [K,F]=assempde(b,p,e,t,c,0,@ntw_pderhs);
    end

    function K=ntw_mat(s)
        K = ntw_assemble(s);
    end

    function F=ntw_rhs(s)
        [K,F] = ntw_assemble(s); %#ok<ASGLU>
    end

    function Av=ntw_matvec(s,v)
        K = ntw_assemble(s);
        Av = K*v;
    end

    function r=ntw_solve(s)
        [K,F] = ntw_assemble(s);
        r = K\F;
    end

P.b = @ntw_rhs;
P.A = @ntw_mat;
P.Av = @ntw_matvec;
P.solve =@ntw_solve;
P.b_const = 0;
P.d = d;
P.N = length(p);
variables(d) = jacobi_parameter;
for i=1:d-1
    variables(i) = jacobi_parameter;
end
P.s = variables;
P.mesh.p = p;
P.mesh.e = e;
P.mesh.t = t;

% this function is 
end

function ft=ntw_pderhs(p,t,u,time)
fn=cos(p(1,:)).*sin(p(2,:));
ft=pdeintrp(p,t,fn');
end


% NTW_COEFF returns a vector with the value of the NTW elliptic coefficient
% given at the cell centers of the mesh.
function ct=ntw_coeff(p,t,s)

d=length(s);
x=p(1,:);
Lc=1/2; Lp=max(1,2*Lc); L=Lc/Lp;

zeta=ones(d+1,1);
zeta(2)=sqrt(sqrt(pi)*L/6);
zeta(3:d+1)=sqrt(sqrt(pi)*L/3)*exp(-0.125*(floor(0.5*(2:d)')*pi*L).^2);

phi=ones(size(x,2),d+1); Lpt=1/Lp;
for i=3:2:d
    phi(:,i)=sin(Lpt*floor(0.5*(i-1))*pi*x');
    phi(:,i+1)=cos(Lpt*floor(0.5*(i-1))*pi*x');
end
if ~mod(d,2)
    phi(:,d+1)=cos(Lpt*floor(0.5*d)*pi*x');
end

v=ones(size(zeta)); v(2:end)=s(:).*zeta(2:end);

% Original NTW coefficient
% cn=0.5+exp(phi*v);

% Constant from coefficient optimized to stay positive for 5 parameters and
% Lc=1/2.
cn=exp(phi*v)-0.85;

ct=pdeintrp(p,t,cn);
end
