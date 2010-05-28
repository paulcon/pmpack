function e = error_estimate(type,X,Y)
%
% e = error_estimate(type,X) computes a particular error estimate based on
% the value of 'type' and an interpolated solution X
% 
%   error_estimate('relerr',X,Y) requires the basis of X to be a subset of
%     the basis of Y.

if nargin<2, error('Not enough inputs.'); end

switch lower(type)
    case 'relerr'
        if ~exist('Y','var') || isempty(Y) 
            error('No reference solution for relative error.'); 
        end
        
        % assumes basis of X is subset of basis of Y
        [Nx,nx]=size(X.coefficients); [Ny,ny]=size(Y.coefficients);
        if nx>ny, T=X; X=Y; Y=T; clear T; t=nx; nx=ny; ny=t; end
        
        e=0;
        for i=1:ny
            match=0;
            for j=1:nx
                if isequal(X.index_set(:,j),Y.index_set(:,i))
                    e=e+norm(X.coefficients(:,j)-Y.coefficients(:,i))^2;
                    match=1;
                    continue
                end
            end
            
            % truncation error
            if ~match 
                e=e+norm(Y.coefficients(:,i))^2;
            end
        end
        e=sqrt(e)/Nx;
    case 'mincoeff'
        e=minimum_coefficient(X);
    case 'resid'
        e=residual_error_estimate(X);
    otherwise
        error('Unrecognized option: %s',type);
end
    