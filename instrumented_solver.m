function S = instrumented_solver(solver,varargin)
%
% R=chol(A(s_mid),'lower');
% pcon1=@(x) galerkin_preconditioner(R',x);
% pcon2=@(x) galerkin_preconditioner(R,x);
% S = instrumented_solver('pcg','tol',1e-8,'maxiter',100,...
%   'precond',{pcon1,pcon2})
% X=spectral_galerkin(Ax,b,s,p_order,'Solver',S.solver);
%
% S.solver - the solver call
% S.stats - results from the solver
% S.histguess - a guess at the time history based on matvec times
%
% Example:
%   Copied from Matlab's PCG, but with increased size
%   n1 = 101; A = gallery('moler',n1);  b1 = A*ones(n1,1);
%   tol = 1e-6;  maxit = 50;  M = diag([(n1-1)/2:-1:1 1 1:(n1-1)/2]);
%   P = chol(M);
%   S = instrumented_solver('pcg','tol',tol,'maxiter',maxit,'precond',{P',P});
%   y=S.solver(A,b1);
%   [flag,resvec] = S.stats()
%   h=S.histguess(); semilogy(h(:,2),h(:,1),'.-'); 
%   xlabel('time'); ylabel('relres');

opts = struct(...
    'tol',1e-8,...
    'precond', [],...
    'restart',15,... % restart for gmres
    'maxiter',500);
for ai=1:2:length(varargin)
    if isfield(opts,varargin{ai}), opts.(varargin{ai})=varargin{ai+1}; end
end

t0 = tic;
Afun = [];
nmult = 0;
matvectimes = zeros(10*opts.maxiter,1);
resvec = [];
relres = Inf;
flag = 1;
histguess = [];
if isempty(opts.precond), opts.precond={[],[]}; end

    function y=Ax(x)
        if isa(Afun,'function_handle')
            y = Afun(x);
        else
            y = Afun*x;
        end
        dt=toc(t0);
        nmult=nmult+1;
        matvectimes(nmult)=dt;
    end

    function x=pcg_solver(A,b)
        Afun = A;
        [x,flag,relres,iter,resvec]=pcg(...
            @Ax,b,opts.tol,opts.maxiter,opts.precond{:});
        matvectimes = matvectimes(1:nmult);
        
        histguess = [resvec/norm(b) matvectimes(1:length(resvec))];
        if length(resvec)+1==length(matvectimes)
            % the final iterate took two-matvecs for checking
            histguess(end,2) = matvectimes(end);
        elseif length(resvec)<length(matvectimes)
            % there were extra checks intersperced
            % TODO do a better interpolation
            histguess(end,2) = matvectimes(end);
        end
    end

    function x=gmres_solver(A,b)
        Afun = A;
        [x,flag,relres,iter,resvec]=gmres(...
            @Ax,b,opts.restart,opts.tol,opts.maxiter,opts.precond{:});
        matvectimes = matvectimes(1:nmult);
       
    end

    function x=minres_solver(A,b)
        Afun = A;
        [x,flag,relres,iter,resvec]=minres(...
            @Ax,b,opts.tol,opts.maxiter,opts.precond{:});
        matvectimes = matvectimes(1:nmult);
        
        histguess = [resvec/norm(b) matvectimes(1:length(resvec))];
        if length(resvec)+1==length(matvectimes)
            % the final iterate took two-matvecs for checking
            histguess(end,2) = matvectimes(end);
        elseif length(resvec)<length(matvectimes)
            % there were extra checks intersperced
            % TODO do a better interpolation
            histguess(end,2) = matvectimes(end);
        end
    end

    function [flag_,resvec_]=solver_stats()
        flag_=flag;
        resvec_=resvec;
    end

    function matvectimes_=timing()
        matvectimes_ = matvectimes; 
    end

    function hist_=hist()
        hist_=histguess;
    end

switch(solver)
    case 'pcg'
        S.solver = @pcg_solver;
        
    case 'gmres'
        S.solver = @gmres_solver;
        
    case 'minres'
        S.solver = @minres_solver;
        
    case 'minres'
        S.solver = @minres_solver;
end

S.stats = @solver_stats;
S.matvectimes = @timing;
S.histguess = @hist;
        

end
