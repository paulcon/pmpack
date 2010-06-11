%% Polynomial expansions
% The key building block of the parametized matrix package (PiMP) is the
% expansion of a vector function in an orthogonal polynomial basis.  This
% demonstration shows how we work with these expansions using our code.
%


%% Single variable polynomial vectors
% A polynomial vector is just like any typical vector, except that the
% elements are polynomial functions.  To get a Matlab vector back from 
% a polynomial vector, we need to evaluate the polynomial at a point.
%
% In PiMP, all polynomial vectors are presented by an expansion in an
% orthogonal polynomial basis.  Let's see how this works in a simple 
% case by defining a polynomial vector.

s = @(x) [2*x + x^2; -x^2 + 2];
% fplot takes vectors the ``wrong'' way, so this function to makes it work
trans2fplot = @(x) s(x)'; 
% plot our vector polynomial 
fplot(trans2fplot,[-1,1]); legend('s1','s2');

%%
% The anonymous function represents $s_1(x)$ and $s_2(x)$ directly with
% their polynomials.  In PiMP, we don't use these representations and
% instead use a basis of orthogonal polynomials.  First, let's see how 
% PiMP would represent the vector polynomial s.

X.x = [legendre_param()];
X.I = [0 1 2];
X.U = [1/3 2/sqrt(3) sqrt(4/45)  ; 2-1/3 0 -sqrt(4/45) ];

% The transpose on the evaluate_expansion is for plotting with fplot
fplot(@(x) evaluate_expansion(X,x)', [-1,1]); legend('r1','r2');

%%
% Although the scale on the plot is slightly different, we hope you see
% these are the same functions.  This methodology may seem needlessly
% complicated for simple functions.  We agree.  Working with the
% representations in orthogonal polynomials, however, gives us the
% mathematical tools we need to write codes that work when we don't have 
% nice closed form solutions.

%% 
% Let's continue by exploring the previous example line-by-line.

%%
% First, we must tell PiMP what kind of parameters we wish to use for our
% function.  The Legendre parameter is the basic parameter... it literally
% means "x" 

%%
% The orthonormal Legendre polynomials over [-1,1] are
%
% $$P_0(x) = 1$$
% $$P_1(x) = \sqrt{3} x $$
% $$P_2(x) = \sqrt{45/4} x^2 - \sqrt{5/4} $$
%


%% 
% One neat side effect is that we can easily plot the orthonormal
% polynomials.

k = 8; % plot the first 10 Legendre polynomials
X.x = [legendre_param()];
X.I = 0:k-1;
X.U = eye(k);
fplot(@(x) evaluate_expansion(X,x)', [-1,1]); legend(cellstr(num2str((X.I)')));

%%
% If we use enough polynomials, it makes a really cool picture.

k = 100; % plot the first 10 Legendre polynomials
X.x = [legendre_param()];
X.I = 1:k;
X.U = eye(k);
[x,y]=fplot(@(x) evaluate_expansion(X,x)', [-1,1]); 
plot(x,y,'b.','MarkerSize',0.5);
axis off; ylim([-3,3]);


%% Multi-variate polynomial vectors
% What happens with multivariate functions?

%% 
% One neat side effect is that we can easily plot the orthonormal
% polynomials.

npts = 50;
k = 3; % plot the first 10 Legendre polynomials
X.x = [legendre_param() legendre_param()];
X.I = [0 0 0 1 1 1 2 2 2; 0 1 2 0 1 2 0 1 2];
X.U = eye(k*k);
[x1,x2] = meshgrid(-1:2/npts:1, -1:2/npts:1);
Zs = arrayfun(@(x) zeros(size(x1)), 1:k*k, 'UniformOutput',false);
for i=1:numel(x1)
    Z = evaluate_expansion(X,[x1(i),x2(i)]);
    for j=1:k*k, 
        Zs{j}(i) = Z(j);
    end
end
figure(1); clf; 
for i=1:k*k
    subplot(k,k,i); mesh(x1,x2,Zs{i}); colormap(bone);
end