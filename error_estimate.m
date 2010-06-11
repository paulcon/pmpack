function e = error_estimate(type,X,Y)
%ERROR_ESTIMATE Estimates the error of a given polynomial approximation.
%
% e = error_estimate(type,X) 
% e = error_estimate('relerr',X,Y) 
%
% The required input 'type' is a string that can take values 'MinCoeff',
% 'RelErr', or 'Resid'. The required input 'X' and optional input 'Y' are 
% structs containing the elements of the polynomial approximation. 
%
% Accepted values of the input 'type' and their associated error
% estimators. 
%   'MinCoeff': Computes the average of the magnitudes of the coefficients
%               associated with the basis polynomials of the two largest
%               degrees. 
%
%   'Resid':    Uses the equation operators A(s) and b(s) to compute a
%               residual error estimate. 
%
%   'RelErr':   Computes the difference between the solution 'X' and a
%               reference solution 'Y'. 
% 
% 
% References:
%   Constantine, P.G., Gleich, D.F., Iaccarino, G. 'Spectral Methods for
%       Parameterized Matrix Equations'. arXiv:0904.2040v1, 2009.
%
% Copyright 2010 David F. Gleich (dfgleic@sandia.gov) and Paul G. 
% Constantine (pconsta@sandia.gov)


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
    