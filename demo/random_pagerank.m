%% Random Alpha PageRank
% See Gleich and Constantine, Random Alpha PageRank, for details about this
% method.  PageRank is a method to evaluate the importance of nodes in a
% directed graph.  On the graph, we follow a random walk where at each
% step,
%   with probability a, we follow a randomly chosen link from the
%      current node
%   with probability (1-a), we jump to a randomly chosen node.
% (Both random choices are uniform over all the possibilities.)
% The stationary distribution of this Markov process is the PageRank
% vector.  We can compute the stationary distribution by solving a linear
% system, 
%
%  $$ (I - aP^T)x = \frac{1-a}{n} e $$
%
% where P is the transition matrix for a random walk on the graph.
% In Random Alpha PageRank, we place a with a random variable, and examine
% the expected value of x.  This computation is solved by a parameterized
% matrix problem.
%
% Usually, we pick a Beta random variable for a -- this corresponds to
% approximating x(a) with a Jacobi parameter.

%% Load a webgraph and setup the linear system
% In this step, we construct the parameterzied linear system that
% correponds to PageRank
load('wb-cs.stanford.mat'); % we load a small web-graph matrix
P = Pcc; % only use the largest strong component -- this fixes a few technical details
A = @(a) (speye(size(P)) - a*P');
b = @(a) (1-a)./size(P,1)*ones(size(P,1),1);
iAb = @(a) A(a)\b(a); % this function solves for a PageRank vector

%% Compute a PageRank vector
% Now, let's compute a single PageRank vector and show the vector
% of ''Ranks'' in sorted order.  This problem corresponds to a surfer
% than will follow a link 85% of the time, or make random jumps
% 15% of the time.  

x85 = iAb(0.85);
semilogy(sort(x85));

%%
% This plot shows that most of the values are small 1e-4 to 1e-3, but
% that a few are large.

%% Setup the parameter
% In random alpha PageRank, we consider a population of surfers where
% the link following probabilities, which are the values of a,
% are generated from a distribution.  
% The parameter s below is a Jacobi parameter to model the randomness in
% the value of a.
% The two shape parameters of the Jacobi parameter give it a slight 
% shift towards larger values of a.

s = [jacobi_parameter(0,1,2,3)];

%%
% There's no good PMPack function to plot the density of such a function.
% The following code draws a beta pdf, which is the density for a Jacobi
% parameter.  
beta_a = 4;  % the code to evaluate the pdf uses beta_a = b+1 = 3+1
beta_b = 3;  %                               and beta_b = a+1 = 2+1
xx = 0:0.01:1;
logkerna = (beta_a-1).*log(xx);
logkernb = (beta_b-1).*log(1-xx);
plot(xx,exp(logkerna+logkernb - betaln(beta_a,beta_b)));

%%
% This figure shows how often surfers will ``draw'' different values of
% a, the link following probabilitiy.

%% Compute the approximation
% At this point, we can compute approximations of the expected value
% of the Random Alpha PageRank function
[X,errz] = pseudospectral(iAb,s,'adapt');

%% Get the expected solution
% The first vector in the coefficients matrix is just the expected
% solution.  It is the expected value of the vector iAb over a drawn from 
% a Beta density.
ex = X.coefficients(:,1);

%%
% Let's compare that against the deterministic solution 
loglog(x85,ex,'.');
xlabel('Deterministic PageRank');
ylabel('Expected PageRank with Jacobi(2,3) parameter');

%%
% The vector has changed quite a bit!  There's a pleasing curve
% in the plot above showing that the values with a Jacobi(2,3)
% are more compressed than their deterministic counterparts.  For
% wild behavior, try a Jacobi(0.2,0.6,0.5,0.5) parameter

%% Plot the first component
% Now we show how a single component of the vector changes as a function
% of a.  
first_component = @(x) x(1);
fplot(@(a) first_component(evaluate_expansion(X,a)), [0,1]);
ylabel('x_1(a)');
xlabel('a');

