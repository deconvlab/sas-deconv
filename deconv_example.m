%% Example - how to use the iPALM code for a CDL problem
clear; clc; 
figure(2); clf;
initpkg;

%% Generate some synthetic data, activation map values are {0,1}.
[A0, X0, Y] = genconvdata(struct('K', 2));

% Mess up the data a little
eta = 5e-2;                             % Add some noise
b0 = 1*rand(size(Y));                   % Add a random constant bias

for i = 1:numel(Y)
    Y{i} = Y{i} + b0(i) + eta * randn(size(Y{i}));
end


%% Set up parameters for iPALM iterations to solve CDL problem. 
% Initial solve
p = size(A0{1});                        % Choose a recovery window size
K = numel(A0);                          % Choose # of atoms to recover
lambda1 = 0.1;                          % Sparsity regularization parameter

xpos = true;                            % Recover a nonnegative activation map
getbias = true;                         % Recover a constant bias

% Reweighting
reweights = 5;                          % number of reweights
lambda2 = 1e-2;                         % lambda for reweighting
eps = 1e-2;                             % reweighting param

% Iterations and updates
maxit = 1e2 * ones(reweights+1,1);      % iterations per reweighting loop
maxit(1) = 2e3;                         % iterations in initial iPALM solve

centerfq = 5e2;                         % frequency to recenter the data
updates = [ 1 10:10:max(maxit) ];       % when to print updates

%% Initialize solver + run some iterations of iPALM

% If K is set to 1, the SBD solver mksbd is more efficient
regularizer = huber(lambda1, xpos);
regularizer = PDRegularizer(regularizer, [], 0.5, 1e-1, [0 Inf]);
solver = cdl_ipalm(Y, p, K, regularizer, getbias);     

update_script = 'cdl_update.m';
figure(1);  subplot(3,2,[1 3]);  
imagesc(Y{1});
title('Original Observation');
figure(2);

% Run iterations
%profile on;
reweight_loop;
%profile off; profile viewer;
