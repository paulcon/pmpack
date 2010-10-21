%% Quadrature and parameters
% PMPack contains all sorts of helpful routines for evaluating and
% computing univariate and multivariate quadrature rules.  
%
% In this demo, we'll see the functions:
% * parameter
% * jacobi_matrix
% * jacobi_eigenvecs
% * gaussian_quadrature
%
% To get started, we need understand where quadrature rules "exist" in
% PMPack.  In brief, they exist with a parameter.

x = parameter();  % consruct a parameter over [-1,1] (this is the default)
[p,w] = gaussian_quadrature(x,5) % build a quadrature rule for that parameter

%%
% *Important Note* 
% *The Quadrature weights in PMPACK are designed to integrate to 1*
% Thus, to get the
% "true" quadrature weights, we must scale the result by two.  To see this,
% consider the integral of the function 1 over the interval [-1,1]
%
% $$\int_{-1}^1 1 \, dx$$

[p,w] =  gaussian_quadrature(parameter,5); % remake the quadrature rule
sum(1.*w)

%%
% But the integral should be two!  That's because we are normalizing the
% weights differently.  To correct, multiply the weight by two.  Or in
% general, multiply by the size of the region of integration.  
sum(1.*w*2)

%% 
% Sticking a quadrature rule into the parameter might seem a bit strange,
% but it's rather handy when we want to combine information from different
% types of parameters, such as when generating multivariate quadrature
% rules.

x = parameter();
y = parameter('legendre',1,2); % construct a parameter over [1,2]
[p,w] = gaussian_quadrature([x,y],3)

%%
% What we get back is the tensor product quadrature rule now.  Because
% all the right information is associated with a parameter, then our codes 
% return the correct points for using a quadrature rule with mixed
% parameters.

%% The parameters themselves
% We now breifly demonstrate some of the details of the parameter objects.
% A parameter itself is a rather simple object

x = parameter()

%%
% It is just a matlab struct with a few different fields.  The name field
% identifies the type.
% In this case, we created a legendre_parameter or a jacobi(0,0) parameter.
% These are equivalent parameters from a quadrature and othogonal
% polynomial point of view.  
% The l and r fields identify the left and right endpoints of the
% parameter.  
% The most interesting field is the "recur" field.  This field is a
% function that returns the recurrence coefficients for the orthogonal
% polynomials associated with this parameter.

ab = x.recur(5)
a = ab(:,1);
b = ab(:,2);

%%
% These recurrence coefficients define the orthogonal polynomials
% associated with a parameter:
%  
% $$P_{-1}(t) = 0$$
%
% $$P_{0}(t) = 1$$
%
% $$\sqrt{b(2)} P_{1}(t) = (t - a(2)) P_{0}(t) - \sqrt{b(1)} P_{-1}(t)$$
%
% $$\sqrt{b(k+2)} P_{k+1}(t) = (t - a(k+1)) P_{0}(t) - \sqrt{b(k+1)}
% P_{k-1}(t)$$
%
% We can encode all of these relationships into a tridiagonal Jacobi
% martrix.  Let $\pi(t)$ be the vector $[P_0(t) \, P_1(t) \, \ldots
% \, P_k(t)]^T$.  Let
%
% $$J = \left[ \begin{array}{ccc} a(1) & \sqrt{b(2)} \\
%                              \sqrt{b(2)} & a(2) & \ddots \\
%                               & \ddots & \ddots \end{array} \right]$$
% 
% Then $t \pi(t) = J pi(t) + \sqrt{b(k+1)} P_{k+1}(t).$
%
% We let you create the matrix J, known as the Jacobi matrix, easily:

J = jacobi_matrix(x,4)

%% 
% Recall the nodes of an n-point quadrature rule is given by the zeros of the
% n-th orthgonal polynomial.  To find these, observe:
%
% $$t \pi(t) = J pi(t) + \sqrt{b(n)} P_{n}(t)$$
%
% where $P_n(t)$ only appears in the final term.  If $P_n(t)$ is 0,
% then 
%
% $$t \pi(t) = J pi(t).$$
%
% This equation is just the eigenvalue eigenvector equation for the matrix
% $J.$  Hence, the eigenvalues of $J$ are the zeros of $P_n(t)$.


%%
% The quadrature points are given by the eigenvalues of this matrix

d = eig(J)

%%
% The quadrature weights are given by the first element of the eigenvectors
% squared.
[V,D] = eig(J)
w = V(1,:).^2

%% 
% We have one more convinence function, which returns the matrix V

V2 = jacobi_eigenvecs(x,4)
norm(V-V2)

%% Convergence
% Let's use the quadrature tools to study the following function

f = @(x) sin(exp(x)).^100;

%%
% First, a plot:

fplot(f,[-1,1])

%%
% To ascertain the accuracy of our quadrature routines, we computed a
% rather accurate answer in Mathematica.
trueval = 0.1598249861160262281333129409016500342546674544621735068837516103348;

%%
% Let us see how we converge to this answer
ks = 1:5:100;
errs = zeros(length(ks),1);
for ki=1:length(ks)
    k = ks(ki);
    [p,w] = gaussian_quadrature(legendre_parameter(),k);
    val = 0;
    for pi=1:k
        val = val + f(p(pi))*2*w(pi); % notice the two here to correct for 
                                      % the weight difference
    end
    errs(ki) = abs(val-trueval);
end
semilogy(ks,errs,'o-');




