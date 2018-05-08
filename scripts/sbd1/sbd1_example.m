%% Example - how to use the iPALM code for a CDL problem
clear; clc; 
run('../../initpkg.m');

%% Generate some synthetic data, activation map values are {0,1}.
% Kernel
p = 20;                         % Size of the (short) kernel

a0 = normpdf(1:p, (p+1)/2, p/10)';
a0 = a0/norm(a0);

% Activation / observation
m = 512;                        % Observation size
theta = 1e-1;                   % Bernoulli (sparsity) coefficient
dist = @(m,n) rand(m,n);        % Distribution of activations

x0 = (rand(m,1) <= theta) .* dist(m,1);

% Mess up the data a little
eta = 0e-2;                    	% Add some noise

y = cconv(a0, x0, m) + eta * randn(m,1);
%stem(x0,'.'); hold on; plot(y-b0); hold off; xlim([1 m]);


%% Set up parameters for iPALM iterations to solve CDL problem. 
lambda = 1e-1;                  % Sparsity regularization parameter
xpos = true;                   % Recover a nonnegative activation map
getbias = false;                % Recover a constant bias

maxit = 2e2;                    % iterations in initial iPALM solve

%% Initialize solver + run some iterations of iPALM
solver = sbd1_ipalm(y, p, lambda, xpos, getbias);

for i = 1:maxit
    solver = iterate(solver);
end

%% Plot results
subplot(311); plot([y cconv(solver.A, solver.X, m)]); xlim([1 m]);
subplot(312); plot([a0 solver.A]); xlim([1 p]);
subplot(313); stem([x0 solver.X], '.'); xlim([1 m]);
maxdotshift(solver.A, a0)

