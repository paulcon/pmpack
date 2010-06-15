
%% error_estimate
    P = pmpack_problem('twobytwo');
    X = spectral_galerkin(P.A,P.b,P.s,2);
    error_estimate('Resid',X)
    
%% evaluate_expansion
    % display a response surface
    P = pmpack_problem('twobytwo');
    X = pseudospectral(P.solve,P.s,5);
    pts = [-1:0.01:1]; fvals_est=zeros(2,length(pts)); fvals_true=fvals_est;
    for i=1:length(pts)
      fvals_true(:,i) = P.solve(pts(i));
      fvals_est(:,i) = evaluate_expansion(X,pts(i));
    end
    plot(pts,fvals_true,'k-',pts,fvals_est,'r-');
    
%% evaluate_ops
    % plot the first 5 Legendre polynomials
    xx = [-1:0.1:1]; % generate a grid
    P = evaluate_ops(parameter,5,xx); % evaluate on xx
    plot(xx,P'); % plot

%% galerkin_preconditioner
    P = pmpack_problem('elliptic1d');
    s0 = midpoint(P.s);
    R = chol(P.A(s0));
    P1 = @(x) galerkin_preconditioner(R',x);
    P2 = @(x) galerkin_preconditioner(R,x);
    X1 = spectral_galerkin(P.A,P.b,P.s,2,...
          'solver',@(A,b) minres(A,b,1e-5,100,P1,P2));
    X2 = spectral_galerkin(P.A,P.b,P.s,2,... % no preconditioning
          'solver',@(A,b) minres(A,b,1e-5,1500));
      
%% gaussian_quadrature
    f = @(x) sin(pi.*x(1))+cos(pi.*x(2));     % integrate the given function
    s = [parameter(); parameter()];
    n = [3; 4];
    [p,w] = gaussian_quadrature(s,n);
    fint = 0;
    for i=1:size(p,1)
        fint = fint + f(p(i,:))*w(i);
    end
    fint
    quad2d(@(X,Y) sin(pi.*X)+cos(pi.*Y),-1,1,-1,1) % compare with matlab quad
    
%% hermite_parameter
    % see that we have variance 1/2
    [p,w] = gaussian_quadrature(s,2);
    std2=0; for i=1:length(p), std2=std2+p(i).^2*w(i); end, std2
    
%% index_set
    % construct a total order basis
    I = index_set('total order',4,2);
    P = pmpack_problem('twobytwo','dim',2);
    X = spectral_galerkin(P.A,P.b,P.s,I);

%% jacobi_eigenvecs

    s = legendre_parameter();
    Q = jacobi_eigenvecs(s,5);

%% jacobi_matrix

    s = legendre_parameter();
    J = jacobi_matrix(s,5);

%% jacobi_parameter

%% 
    
%% pseudospectral
    A = @(t) [2 t; t 1];                    % 2x2 parameterized matrix
    b = @(t) [2; 1];                        % constant right hand side
    iAb = @(t) A(t)\b(t);
    s = parameter();                        % parameter defined on [-1,1]
    pOrder = 13;                            % degree 13 approximation
    X = pseudospectral(iAb,s,pOrder);
    

%% spectral_galerkin
    A = @(t) [2 t; t 1];                    % 2x2 parameterized matrix
    b = @(t) [2; 1];                        % constant right hand side
    s = parameter();                        % parameter defined on [-1,1]
    pOrder = 13;                            % degree 13 approximation
    X = spectral_galerkin(A,b,s,pOrder);
    
