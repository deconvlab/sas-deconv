%% Example - how to use the iPALM code for a CDL problem
clear; clc; 
run('../../initpkg');

%% Generate some synthetic data, activation map values are {0,1}.
[A0, X0, Y] = genconvdata(struct('K', 1, 'm', 512));

% Mess up the data a little
eta = 0e-2;                             % Add some noise
b0 = 1*rand;                   % Add a random constant bias
Y = Y + b0;

%% Set up parameters for iPALM iterations to solve CDL problem. 
% Initial solve
p = size(A0);                        % Choose a recovery window size
lambda = 0.1;                          % Sparsity regularization parameter

xpos = true;                            % Recover a nonnegative activation map
getbias = true;                         % Recover a constant bias

maxit = 1e4;                         % iterations in initial iPALM solve
update = 10;


%% Initialize solver + run some iterations of iPALM
solver = mksbd(Y, p, lambda, true, true);
costfun = get(solver,'H');
costfun = @costfun.value;

solver = admm_sbd(Y, p, huber(lambda, xpos));
solver.X = ones(size(Y));
solver.rhoA = [20 500 2];  
solver.rhoX = [20 500 1.2];

costs = NaN(maxit, 4);
for i = 1:maxit
    solver = iterate(solver);
    costs(i,1:4) = solver.cost;
    costs(i,5) = costfun({solver.A}, {solver.X}, {solver.b}, []);
    
    if mod(i, update) == 0 
        figure(1);
        
        subplot(231); imagesc(Y);
        title(['Iteration ' num2str(solver.it)]);
        subplot(234); imagesc(cconvfft2(solver.A, solver.X));
        title(['Cost ' num2str(solver.cost(1))]);
        
        subplot(232); imagesc(A0);
        subplot(235); imagesc(solver.A);

        subplot(233); imagesc(X0);
        subplot(236); imagesc(solver.X);
        
        drawnow;
    end
end
solver %#ok<NOPTS>
