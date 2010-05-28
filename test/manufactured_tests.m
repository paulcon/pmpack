%% Manufactured solutions

msgid = 'pmpack:manufacturedTests';

%% one parameter tests (A+sB)*(x1+s*x2) = A*x1 + s*A*x2 + s*B*x1 + s^2 B*x2

%% 
% pseudospectral
for i=1:50
    A = randn(10);
    B = randn(10);
    A=A'*A;
    B=B'*B;
    x1 = randn(10,1);
    x2 = randn(10,1);
    
    Afun = @(s) A + s*B;
    bfun = @(s) (A+s*B)*(x1+s*x2);
    xtrue = @(s) x1+s*x2;
    iAb = @(s) Afun(s)\bfun(s);
        
    ab=floor(16*rand(2,1)); lr=sort(rand(2,1));
    s = jacobi_parameter(lr(1),lr(2),ab(1),ab(2));
    
    X = pseudospectral(iAb,s,2);
    for si=linspace(lr(1),lr(2),15);
        xsi = evaluate_expansion(X,si);
        if norm(xtrue(si)-xsi,'inf')>1e-12
            error(msgid,'pseudo-spectral failed manufactured solution:');
        end
    end
    
    X = spectral_galerkin(Afun,bfun,s,2);
    for si=linspace(lr(1),lr(2),15);
        xsi = evaluate_expansion(X,si);
        if norm(xtrue(si)-xsi,'inf')>1e-12
            error(msgid,'galerkin failed manufactured solution:');
        end
    end
end
            

