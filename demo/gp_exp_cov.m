function C = gp_exp_cov(x1,x2,c,sigma)
% GP_EXP_COV is the two point exponential covariance function with
% correlation length given by the d-vector 'c'. The scalar 'sigma'
% amplifies the covariance.
%
% C = gp_exp_cov(x1,x2,c,sigma)
%
% Both 'x1' and 'x2' must be row vectors of length d. 

if max(size(c))==1, c=c*ones(size(x1)); end
r=(x1-x2)*diag(1./c.^2)*(x1-x2)';
C=sigma*exp(-0.5*r);