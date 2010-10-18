% Parameterized Matrix Package
% Version 1.0 
%
% Approximation methods
% ---------------------
%   pseudospectral          - Pseudospectral approximation of solution to A(s)x(s)=b(s)
%   spectral_galerkin       - Galerkin approximation of solution to A(s)x(s)=b(s)
%
% Utility functions
% -----------------
%   error_estimate          - Estimates the error of a given polynomial approximation
%   evaluate_expansion      - Evaluates a polynomial expansion at a point
%   evaluate_ops            - Evaluate the n-vector of orthonormal polynomials at a point
%   galerkin_preconditioner - 
%   gaussian_quadrature     - Compute a Gaussian quadrature rule
%   index_set               - Builds the multi-indicies for orthogonal polynomials
%   midpoint                - Computes the midpoint of the domain defined by 's'
%   minimum_coefficient     - Type of error estimate for given 'X'
%   residual_error_estimate - Computes residual error estimate
%   sort_bases              - Sorts the multivariate bases by total degree
%   pmpack_problem          - Load a problem structure for a parameterized matrix problem
%
% Parameter functions
% -------------------
%   hermite_parameter       - Construct a parameter with a Gaussian weight function
%   jacobi_eigenvecs        - Construct the eigenvectors of the Jacobi matrix
%   jacobi_matrix           - Construct a Jacobi matrix for a parameter
%   jacobi_parameter        - Construct a parameter with a Jacobi weight function
%   legendre_parameter      - Construct a parameter with a uniform weight function
%   parameter               - Construct a parameter with a uniform weight function
%
% References
% ----------
%   Constantine, P.G., Gleich, D.F., Iaccarino, G. 'Spectral Methods for
%       Parameterized Matrix Equations'. SIMAX, 2010.
%       http://dx.doi.org/10.1137/090755965
%
%   Constantine, P.G., Gleich, D.F., Iaccarino, G. 'A factorization of 
%       the spectral Galerkin system for parameterized matrix equations: 
%       Derivation and applications'. 2010
%
% Copyright 2010 David F. Gleich (dfgleic@sandia.gov) and Paul G. 
% Constantine (pconsta@sandia.gov).

