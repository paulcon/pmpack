function s = hermite_parameter()
% HERMITE_PARAM Construct a parameter with a Hermite weight function
%
% x = hermite_param() generates a new Hermite parameter over [-Inf,Inf] with 
% a Gaussian weight.  This parameter is the canonical integration parameter.
%
% Copyright, Stanford University, 2009
% Paul G. Constantine, David F. Gleich

s.name = sprintf('hermite');
s.recur = @(n) hermite_recur(n);

end

function ab=hermite_recur(n)
% Compute the recurrence coefficients for the Hermite polynomials using
% Gautschi's OPQ software.
ab=r_hermite(n);
ab(1,2)=1;
end