function test_main

galerkin_testscript;
pseudospec_testscript;
utils_testscript;
examples_test;

% more exhaustive tests
pagerank_tests
quadrature_tests
manufactured_tests
problems_tests

if license('test','PDE_Toolbox')
    if license('checkout','PDE_Toolbox')
        test_main_pdetoolbox
    else
        fprintf('Cannot checkout the PDE toolbox for tests\n');
    end
else
    fprintf('Skipping tests that require the PDE toolbox\n');
end
