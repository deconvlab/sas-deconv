clc; clear; %#ok<*NOPTS>
run('../../initpkg.m');
addpath('colormapline');
addpath('../admm');

%% Problem setup
clf; hold off; drawnow;

p = 32;
a0 = randn(p,1);  a0 = a0/norm(a0);

m = 1024;                       % Observation size
theta = 2e-1;                   % Bernoulli (sparsity) coefficient
dist = @(m,n) randn(m,n);       % Distribution of activations
x0 = (rand(m,1) <= theta) .* dist(m,1);

y = cconv(a0, x0, m);

lambda = 0.1;


% Plot phi over the simplex
s = 4;  s = -ceil(p/s):ceil(p/s);
s = s(randperm(numel(s), 3));
s = [0 s(s~=0)];  s = s(1:3);

nsamps = 1e2;
samples = linspace(-1,1.5, nsamps);

a1 = [a0; zeros(m-p,1)];
a2 = circshift(a1, s(3));  a2 = a2(1:p);  a2 = a2/norm(a2);
a1 = circshift(a1, s(2));  a1 = a1(1:p);  a1 = a1/norm(a1);
costs = repmat({NaN(1,nsamps)}, [nsamps 1]);
parfor u = 1:nsamps
    phi = phi_fista(huber(lambda));
    for v = 1:nsamps
        a = a0 + samples(u)*(a1-a0) + samples(v)*(a2-a0); %#ok<PFBNS>
        a = a/norm(a);
        [phi, costs{u}(v)] = evaluate(phi, y, a);
    end
end
costs = cell2mat(costs);

surf(samples, samples, costs);  hold on;  colormap parula
shading interp;  view(2);  drawnow;
xlabel('u');  ylabel('v');

contour(samples, samples, costs);  colormap parula
plot([0 1; 0 0; 0 1]', [0 0; 0 1; 1 0]', 'k', 'LineWidth', 1.2);
plot(1/3, 1/3, 'o');

%% Plot trajectories over hemisphere
maxit = 300;

% Fix initial point
%a = [0.46 -0.68];
%a = [-.69 .054];
a = [.55 .5];
a(3) = sqrt(max(1-sum(a.^2),0));

solvers = {
    sbd1_ipalm(y, 3, lambda, false, false, a) ...
    sbd1_ipalm(y, 3, lambda, false, false, a) ...
    admm_sbd(y(:), [3 1], huber(lambda), a(:)) };

solvers{1}.iterator.alph = 0;
solvers{2}.iterator.alph = 0.9;
solvers{3}.rhoA = [1 30 2];
solvers{3}.rhoX = [1 30 2];
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