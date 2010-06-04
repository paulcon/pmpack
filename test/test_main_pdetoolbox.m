function test_main_pdetoolbox
% Tests that require the pdetoolbox



addpath('../problems')
try
    %% NTW
    P = nobile_tempone_webster();
    [X,r] = spectral_galerkin(P.Av,P.b,P.s,'adapt');
    scatter(P.mesh.p(1,:),P.mesh.p(2,:),1,X.coefficient(:,1));
    
    rmpath('../problems');
catch me
    rmpath('../problems');
    throw(me);
end
