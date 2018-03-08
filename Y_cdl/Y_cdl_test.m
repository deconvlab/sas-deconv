clear; clc; 

dirs = {'../', '../helpers'};
for i = 1:numel(dirs);  addpath(dirs{i});  end

%% Testing SBD reweighting with Y.mat
% Load data and choose p0, K
load('Y_x00_71K_256_blueSpots.mat')
Y = Y/max(abs(Y(:)));
resz = 2;                               % resizing factor
Y = {imresize(Y,resz)};

p = [10 10] * resz;                     % window size
K = 2;                                  % # of kernels to recover

%% Set up parameters for solving the CDL problem. 
% Initial solve
xpos = true;        % recover a nonnegative activation map
getbias = true;     % recover a constant bias
lambda1 = 5e-2;     % sparsity regularization parameter

% Reweighting
rwgts = 5;                              % number of reweights
lambda2 = 5e-3;                         % lambda for reweighting
eps = 1e-2;                             % reweighting parameter

% Maximum iterations
maxit = [1e3; 5e1*ones(rwgts,1)];       % iterations per outer loop                        

% Centering and updates
centerfq = 0e2;                         % how often to recenter the data
updates = [1 10:10:max(maxit)];         % when to print updates

% Initialize solver
[sol, synthesize] = mkcdl(Y, p, K, lambda1, true, true);
%sol.info.debug = true;

%% Run some iterations of iPALM
sols = cell(rwgts+1,1);
costs = cell(rwgts+1,1);  
stime = tic;  profile on;
figure(1);  subplot(3,2,[1 3]);  imagesc(Y{1});
for r = 1:rwgts+1
    if r > 1
        sol = reweight(sol, lambda2, eps);
    else
        figure(2); clf;
    end
    costs{r} = NaN(maxit(r),1);
    
    for i = 1:maxit(r)
        sol = iterate(sol);
        costs{r}(i) = sol.cost;

        % TODO: Centering needs to be fixed
        if centerfq > 0 && mod(i, centerfq) == 0
            [A, X, ~, A_, X_] = center(...
                sol.A{:}, sol.X{:}, ...
                [], sol.A_{:}, sol.X_{:});
            sol.A = {A};  sol.X = {X};
            sol.A_ = {A_};  sol.X_ = {X_};
        end
        
        % Updates and figures
        if ismember(i, updates)
            tmp = ['Iteration ' num2str(sol.it)];
            
            figure(1);
            subplot(3,2,[2 4]); 
                imagesc(synthesize(sol.A, sol.X, sol.b));
            subplot(313); 
            for k = 1:r
                c = costs{k};
                if k == 1;  c = c(50:end);  end
                plot(linspace(0,1,numel(c)), c-min(c));  hold on;
            end
            hold off; title(tmp); drawnow;
            
            figure(2);
            for k = 1:K
                subplot(2,2*K,2*K*(r>1) + 2*k-1); 
                imagesc(sol.A{k});
                subplot(2,2*K,2*K*(r>1) + 2*k);
                imagesc(sol.X{k});
            end
            subplot(2,2*K,2*K*(r>1)+1); title(tmp);
            drawnow;
            
            fprintf(['Iter. %d:%d.  '...
                'Cost %.4e. Elapsed time %.2fs.\n'], ...
                r-1, sol.it, sol.cost, toc(stime));
        end
    end
    disp(' ');
    sols{r} = copy(sol);
end
profile off; profile viewer;

%% Keep only relevant variables
clearvars ans r i k c sol stime tmp updates dirs

solutions = cell(2,1);
solutions{1} = get(sols{1}, 'A', 'X', 'b');
solutions{2} = get(sols{end}, 'A', 'X', 'b');