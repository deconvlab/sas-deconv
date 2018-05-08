%% SBD1 Phase transition
clear; clc; 
run('../../initpkg.m');

%% Generate some synthetic data, activation map values are {0,1}.
% Activation / observation
m = 2048;                       % Observation size
dist = @(m,n) randn(m,n);       % Distribution of activations

% Trials
trials = 10;
p = 50;
tfracs = linspace(0,3,20);
thetas = 10.^tfracs/p;

%% Set up parameters for iPALM iterations to solve CDL problem. 
lambda = 1e-1;                  % Sparsity regularization parameter
xpos = false;                   % Recover a nonnegative activation map
getbias = false;                % Recover a constant bias

%% Initialize solver + run some iterations of iPALM
costs = NaN(numel(thetas), trials);
its = NaN(numel(thetas), trials);
a0 = randn(p,1);  a0 = a0/norm(a0);

for i = 1:numel(thetas)
    theta = thetas(i);
parfor trial = 1:trials
    %a0 = randn(p,1);  a0 = a0/norm(a0);
    xgood = false;
    while ~xgood
        x0 = (rand(m,1) <= theta) .* dist(m,1);
        xgood = sum(x0~=0) >= 1;
    end
    y = cconv(a0, x0, m);
    solver = sbd1_ipalm(y, 3*p-2, lambda, xpos, getbias);

    solver = loop_phasetrans(solver);
    
    costs(i, trial) = maxdotshift(a0, solver.a);
    its(i, trial) = solver.it;
end
end

%%
plot(tfracs, mean(costs,2));
