%% Polynomial expansions
% The key building block of the parametized matrix package (pmpack) is the
% expansion of a vector function in an orthogonal polynomial basis.  This
% demonstration shows how we work with these expansions using our code.
%
% In particular, it demonstrates the codes
%
% * evaluate_expansion
% * evaluate_ops
% * pseudospectral
% 

%% Single variable polynomial vectors
% A polynomial vector is just like any typical vector, except that the
% elements are polynomial functions.  To get a Matlab vector back from 
% a polynomial vector, we need to evaluate the polynomial at a point.

s = @(x) [2*x + x^2; ...
          -x^2 + 2];
  
%%
s(0)
%%      
s(1)
%%
s(2)

%%
% We can also plot the polynomials themselves.

% fplot takes vectors the ``wrong'' way, so this function to makes it work
trans2fplot = @(x) s(x)'; 
% plot our vector polynomial 
fplot(trans2fplot,[-1,1]); legend('s1','s2');

%%
% The anonymous function represents $s_1(x)$ and $s_2(x)$ directly with
% their polynomials.  In PMPack, we don't use these representations and
% instead use a basis of orthogonal polynomials.  First, let's see how 
% PMPack would represent the vector polynomial s 

X.variables = [legendre_parameter()];
X.index_set = [0 1 2];
X.coefficients = [1/3 2/sqrt(3) sqrt(4/45)  ; 2-1/3 0 -sqrt(4/45) ];

%%
% This new variable X is the PMPack representation of the polynomial vector
% s.

evaluate_expansion(X,0)
%%
evaluate_expansion(X,1)
%%
evaluate_expansion(X,2)

%%
% And now let's draw the polynomials themselves.

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
% First, we must tell PMPack what kind of parameters we wish to use for our
% function.  The Legendre parameter is the basic parameter... it literally
% means "x" when your function is defined over [-1,1].

X.variables = [legendre_parameter()];

%%
% The orthonormal Legendre polynomials over [-1,1] are
%
% $$P_0(x) = 1$$
%
% $$P_1(x) = \sqrt{3} x $$
%
% $$P_2(x) = \sqrt{45/4} x^2 - \sqrt{5/4} $$
%
% To plot the polynomials $s_1(x)$ and $s_2(x)$, we need to figure out
% how they are expressed in these as multiples of these three polynomials.
% We only need three polynomials, and thus, we create the index set:

X.index_set = [0 1 2];

%%
% This variable tells pmpack that we want to work with the 0th, 1st, and
% 2nd orthonormal Legendre polynomials.
%
% We'll leave the next step to you to check, but
%
%  $$s_1(x) = 2*x + x^2 = 1/3*P_0(x) + 2/sqrt(3)*P_1(x) + sqrt(4/45)*P_2(x)$$
%
% and
%
%  $$s_2(x) = -x^2 + 2 = 5/3*P_0(x) - sqrt(4/45)*P_2(x)$$
%
% in which case, we can tell pmpack that we want these combinations
% of P_0, P_1, and P_2 with

X.coefficients = [1/3 2/sqrt(3) sqrt(4/45); % for s_1
                  2-1/3 0 -sqrt(4/45) ];    % for s_2
              
%%
% The relationship between X.coefficients and X.index_set is:
%   x(1) = sum over i
%            X.coefficients(1,i)*orthogonal_polynomial(X.index_set(i))
% so the index set controls which orthogonal polynomial, and the
% the coefficients control how much of that polynomial.
   
%% 
% The command
%   evaluate_expansion
% takes in a representation of a function in orthogonal polynomials 
% and returns the value of that function at a point.
% In this case,

evaluate_expansion(X,0.5)

% is doing the same thing as

s(0.5)

%% 
% The last thing we did in the previous example was draw a plot.  This
% just used fplot to call evaluate_expansion at many point:

fplot(@(x) evaluate_expansion(X,x)', [-1,1]); legend('r1','r2');


%% Drawing orthogonal polynomials
% One neat feature is that we can easily plot these orthonormal
% polynomials.  This occurs because we write each of our functions
% in a basis of orthogonal polynomials.  To draw the polynomials
% themselves, then, we just use the _identity_ matrix as the combination
% of polynomials.  Let's see an example

k = 8; % plot the first 8 Legendre polynomials
X.variables = [legendre_parameter()];
X.index_set = 0:k-1;  
X.coefficients = eye(k); % the matrix U tells pmpack how much weight to put on
              % each orthogonal polynomial, in this case, we are just
              % weighting the diagonal, which will draw all the
              % polynomials.

              
fplot(@(x) evaluate_expansion(X,x)', [-1,1]); legend(cellstr(num2str((X.index_set)')));

%%
% If we use enough polynomials, it makes a really cool picture.

k = 100; % plot the first 10 Legendre polynomials
X.variables = [legendre_parameter()];
X.index_set = 1:k;
X.coefficients = eye(k);
[x,y]=fplot(@(x) evaluate_expansion(X,x)', [-1,1]); 



%%
% But there is an easier function for this task.  The 
% code evaluation_ops does exactly the same thing!

fplot(@(x) evaluate_ops(legendre_parameter(),8,x)', [-1,1]);

%%
% The Jacobi polynomials
fplot(@(x) evaluate_ops(jacobi_parameter(-1,1,2,3),8,x)', [-1,1]);
ylim([-5,5]);

%%
% The Hermite polynomials
fplot(@(x) evaluate_ops(hermite_parameter,8,x)', [-2,2]);
ylim([-5,5]);


%% Multi-variate polynomial vectors
% What happens with multivariate functions?
% In the previous demostrations, everything used univariate polynomials.
% It's easy enough in Matlab to create a multi-variable polynomial vector.

% These are defined elementwise, then grouped so we can do plots below.
s = @(x1,x2) [x1.*x2 + 1; x1.^2 + x2.^2];
%%
% These are a slightly tricker to plot, but not too bad.  This is mainly
% due to Matlab's idiocyncracies.  
vecplot = @(X,Y,f) mesh(X,Y,reshape(f(X(:)',Y(:)'),[size(X,1),size(X,2)]));
[X,Y] = meshgrid(-1:0.1:1,-1:0.1:1);
ith_component = @(x,i) x(i,:);
subplot(1,2,1); vecplot(X,Y,@(x,y) ith_component(s(x,y),1));
subplot(1,2,2); vecplot(X,Y,@(x,y) ith_component(s(x,y),2));

%% 
% Rather that figuring out the coefficients of the expansion ourselves,
% we can use the pseudo-spectral code to do it!
% We won't get into why exactly this works here, but this is equivalent to
% the problem that the pseudospectral code solves.  

X = pseudospectral(@(t) s(t(1),t(2)), [parameter(),parameter()], 2);

%%
% The matrix X now contains our expansion in multi-variate orthogonal
% polynomials
X

%%
% X.variables is always just a copy of [parameter(),parameter()].

%%
% Let's look at these results
X.index_set

%%
% For univariate functions, index_set is just a list of orthogonal
% polynomials.  For multivariate functions, the orthogonal polynomials are
% tensor products of 1d-orthgonal polynomials.  Thus,
% X.index_set(:,1) = [0,0] 
% means to use P_0(x1)*P_0(x2) 
% and
% X.index_set(:,2) = [0,1]
% means P_0(x1)*P_1(x2)
% These index sets get complicated rather quickly, but we need them to
% express the polynomials.

%%
% Let's look at the coefficients
X.coefficients

%%
% So what this says is that 
% s1(x1,x2) = 1*P_0(x1)*P_0(x2) + 1/3*P_1(x1)*P_1(x2)
% s2(x1,x2) = 2/3*P_0(x1)*P_0(x2) + sqrt(4/45)*P_0(x1)*P_2(x2)*P
%               + sqrt(4/45)*P_2(x1)*P_0(x2)

%%
% Let's see what this plot looks like!  This plot is a little more annoying
% because evaluate expansion only works at a single point right now.
[GX,GY] = meshgrid(-1:0.1:1,-1:0.1:1);
% gather data
Z1 = zeros(length(GX));
Z2 = Z1;
for i=1:numel(GX)
    z = evaluate_expansion(X,[GX(i),GY(i)]);
    Z1(i) = z(1);
    Z2(i) = z(2);
end
subplot(1,2,1); mesh(GX,GY,Z1);
subplot(1,2,2); mesh(GX,GY,Z2);
    
              
%% Multivariate orthogonal polynomials
% Again, we can plot the orthogonal polynomials too!

npts = 50;
k = 3; % plot the first 10 Legendre polynomials
X.variables = [legendre_parameter() legendre_parameter()];
X.index_set = [0 0 0 1 1 1 2 2 2; 0 1 2 0 1 2 0 1 2];
X.coefficients = eye(k*k);
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