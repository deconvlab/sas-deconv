clear; clc; 
figure(2); clf;

dirs = {'../', '../helpers'};
for i = 1:numel(dirs);  addpath(dirs{i});  end

%% Testing SBD reweighting with Y.mat
load('Y.mat');          % Load data and choose p0
%load('/home/yenson/Desktop/Y_x015_10meV_6K_128_singleBoundary.mat')
Y = Y/max(abs(Y(:)));
resz = 2;
Y = imresize(Y,resz);
p = [10 10] * resz;

%% Initialize an iPALM iterator for the SBD problem. 
xpos = true;        % Recover a nonnegative activation map
getbias = true;     % Recover a constant bias
lambda1 = 1e-1;     % Sparsity regularization parameter

solver = mkcdl(Y, p, 1, lambda1, true, true);     

%% Run some iterations of iPALM
reweights = 5;                          % number of reweights
lambda2 = 1e-2;                         % lambda for reweighting
eps = 1e-2;                             % reweighting param

maxit = 1e2 * ones(reweights+1,1);      % iterations per reweighting loop
maxit(1) = 2e3;                         % iterations in initial iPALM solve

centerfq = 0e2;                         % frequency to recenter the data
updates = [ 1 2:2:10 ...                % when to print updates
            50:50:200 ...
            400:200:max(maxit)];

solvers = cell(reweights+1,1);
costs = cell(reweights+1,1);  
stime = tic;  %profile on;
for r = 1:reweights+1
    if r > 1
        solver = reweight(solver, lambda2, eps);
    end
    costs{r} = NaN(maxit(r),1);
    
    for i = 1:maxit(r)
        solver = iterate(solver);
        costs{r}(i) = solver.cost;

        if centerfq > 0 && mod(i, centerfq) == 0
            [A, X, ~, A_, X_] = center(...
                solver.A{:}, solver.X{:}, ...
                [], solver.A_{:}, solver.X_{:});
            solver.A = {A};  solver.X = {X};
            solver.A_ = {A_};  solver.X_ = {X_};
        end
        
        if ismember(i, updates)
            figure(1); imgupdate(Y, solver);  drawnow;
            
            figure(2);             
            if r < 2
                subplot(221); imagesc(abs(solver.A{:}));
                subplot(223); imagesc(solver.X{:});
            else
                subplot(222); imagesc(abs(solver.A{:}));
                subplot(224); imagesc(solver.X{:});
            end
            
            fprintf(['Iter. %d:%d.  '...
                'Cost %.4e. Elapsed time %.2fs.\n'], ...
                r-1, solver.it, solver.cost, toc(stime));
        end
    end
    disp(' ');
    solvers{r} = copy(solver);
end
%profile off; profile viewer;