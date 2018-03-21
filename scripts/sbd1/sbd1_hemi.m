clc; clear; %#ok<*NOPTS>
clf; hold off; drawnow;
run('../../initpkg.m');
addpath('colormapline');
addpath('../admm');

%% Problem setup
a0 = [1 1 0.2]/sqrt(3);

m = 1024;                        % Observation size
theta = 2e-1;                   % Bernoulli (sparsity) coefficient
dist = @(m,n) randn(m,n);       % Distribution of activations
x0 = (rand(m,1) <= theta) .* dist(m,1);

y = cconv(a0, x0, m);

lambda = 0.1;

%% Plot phi over the hemisphere
samples = sample_hemi3(40,200);

costs = repmat({NaN(1,size(samples,2))}, [size(samples,1) 1]);
parfor i = 1:size(samples,1)
    phi = phi_fista(huber(lambda));
    for a = 1:size(samples,2)
        [phi, costs{i}(a)] = evaluate(phi, y, samples(i,a,:));
    end
end
costs = cell2mat(costs);

surf(samples(:,:,1), samples(:,:,2), costs);  hold on;  colormap jet
contour(samples(:,:,1), samples(:,:,2), costs);  colormap jet
shading interp;  view(2);
xlabel('x');  ylabel('y');

%% Plot trajectories over hemisphere
maxit = 300;

% Fix initial point
a = [0.46 -0.68];
%a = [-.69 .054];
a = [a sqrt(1-sum(a.^2))];

solvers = {
    sbd1_ipalm(y, 3, lambda, false, false, a) ...
    sbd1_ipalm(y, 3, lambda, false, false, a) ...
    admm_sbd(y(:), [3 1], huber(lambda), a(:)) };

solvers{1}.iterator.alph = 0;
solvers{2}.iterator.alph = 0.9;
solvers{3}.rhoA = [1 20 2];
solvers{3}.rhoX = [1 20 2];
solvers{3}.rhoZ = [1 0 1]

phi = arrayfun(@(~) evaluate(phi_fista(huber(lambda)), y, a), ...
    1:numel(solvers), 'UniformOutput', false);

samples = repmat({[a; NaN(maxit,3)]}, [1 numel(solvers)]);
costs = repmat({[phi{1}.costs(end); NaN(maxit,1)]}, [1 numel(solvers)]);
for i = 2:maxit+1
    for j = 1:numel(solvers)
        solvers{j} = iterate(solvers{j});
        samples{j}(i,:) = solvers{j}.A;
        [phi{j}, costs{j}(i)] = evaluate(phi{j}, y, solvers{j}.A);
    end
end

% Plotting
lgd = {'PALM', 'iPALM', 'ADMM'};
colors = [0 1 0; 1 0 1; 1 .5 .3];
sym = {'x', 'o', '^'};

h = cell(1,2);
for j = 1:numel(solvers)
    h{1} = [h{1}, colormapline(...
        samples{j}(:,1), samples{j}(:,2), costs{j}+1, ...
        min(log10(linspace(1,10,maxit)),1)'*colors(j,:) ...
    )];

    colormapline(samples{j}(:,1), samples{j}(:,2), [], ...
        min(log10(linspace(1,10,maxit)),1)'*colors(j,:));
end

i = 1:10:maxit+1;
for j = 1:numel(solvers)
    h{2} = [h{2}, plot3(...
        samples{j}(i,1), samples{j}(i,2), costs{j}(i)+1, ...
        sym{j}, 'LineWidth', 1, 'Color', colors(j,:) ...
    )];

    plot(samples{j}(i,1), samples{j}(i,2), ...
        sym{j}, 'LineWidth', 1, 'Color', colors(j,:));
end

legend(h{2}, lgd(1:numel(solvers)));