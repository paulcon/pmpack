function P = elliptic_func(varargin)
% ELLIPTIC_FUNC Construct a parameterized elliptic PDE.
%
% This function requires the PDE toolbox.  It also produces a warning when
% we generate the mesh, you can safely ignore this warning.
%
% P = elliptic_func(); 
% P = elliptic_func(...); 
%
% Constructs a problem struct with the parameterized matrix interfaces for
% an elliptic PDE in two spatial dimensions discretized on an irregular
% mesh with parameterized, spatially dependent coefficients representing a
% truncate Karhunen-Loeve expansion of a random field.
%
% Outputs:
% The output P is a struct with the following fields:
%   P.A:     A parameterized matrix interface.
%   P.Av:    A parameterized matrix-vector interface.
%   P.b:     A parameterized right hand side.
%   P.N:     The size of the matrix.
%   P.d:     The number of parameters.
%   P.solve: A parameterized solution.
%   P.s:     The vector of parameters for this problem.
%   P.mesh:  A struct with mesh details (points, edges, triangles).
%   P.KL:    A struct containing the eigenvectors and singular values of
%            the KL expansion.
%
% Optional inputs:
% To specify optional inputs, use the 'key/value' format. For example, to
% set the correlation length 'corr' to 0.5, include 'corr',0.5 in the 
% argument list. See the examples below for more details.
%   corr:   Correlation length of the random field model of the elliptic
%           coefficients (Default 1)
%   scale:  A scaling parameter for the magnitude of the random field.
%           (Default 1)
%   trunc:  Truncation level of the random field. determines the dimension
%           of the input parameter space. DON'T MAKE THIS TOO LARGE.
%           (Default 4)
%   refine: Set this to 1 to refine the spatial mesh. Unrefined mesh has
%           177 nodes. The refined mesh has 665 nodes. (Default 0)
%   
%
% Example:
%   P = elliptic_func('trunc',3);
%   [X,r] = spectral_galerkin(P.Av,P.b,P.s,2);
%   pdesurf(P.mesh.p,P.mesh.t,X.coefficient(:,1)); 
%   view(2); colorbar;
%

% Copyright, Stanford University, 2009-2010
% Paul G. Constantine, David F. Gleich

corr = 1;
scale = 1;
d = 3;
refine = 0;

for i=1:2:nargin
    switch lower(varargin{i})
        case 'corr'
            corr = varargin{i+1};
        case 'scale'
            scale = varargin{i+1};
        case 'trunc'
            d = varargin{i+1};
        case 'refine'
            refine = varargin{i+1};
        otherwise
            error('Unrecognized option: %s',varargin{i});
    end
end


% Geometry and boundary conditions
[p,e,t] = initmesh('squareg');
if refine, [p,e,t] = refinemesh('squareg',p,e,t); end
fprintf('Created a mesh with %i nodes.\n', length(p));
b = 'squareb1';

% Random field model of the elliptic coefficients
cv = @(x1,x2) gp_exp_cov(x1,x2,corr,scale);
[F,KL] = randomfield(cv,p','trunc',d);
U = KL.bases*diag(KL.sv);

% Right hand side
f_node=(cos(p(1,:)).*sin(p(2,:)))';
f_tri=pdeintrp(p,t,f_node);

% set these as nested functions so they share all the loaded data
    function [K,F] = assemble(s) 

        c_node = exp(U*s(:));
        c_tri = pdeintrp(p,t,c_node);

        % Solve the system.
        [K,F]=assempde(b,p,e,t,c_tri,0,f_tri);
    end

    function A = matfun(s)
        A = assemble(s);
    end

    function Av = matvecfun(s,v)
        K = matfun(s);
        Av = K*v;
    end

    function r = solve(s)
        [K,F] = assemble(s);
        r = K\F;
    end

[~,F] = assemble(zeros(1,d));
P.b = @(s) F; % parameter independent forcing function
P.A = @matfun;
P.Av = @matvecfun;
P.solve = @solve;
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
P.KL = KL;

end
