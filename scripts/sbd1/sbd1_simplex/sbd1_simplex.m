clc; clear; %#ok<*NOPTS>
run('../../../initpkg.m');
addpath('colormapline');
addpath('../../admm');
addpath('../');

%% Problem setup
p = 32;
a0 = randn(p,1);  a0 = a0/norm(a0);

m = 1024;                       % Observation size
theta = 2e-1;                   % Bernoulli (sparsity) coefficient
dist = @(m,n) randn(m,n);       % Distribution of activations
x0 = (rand(m,1) <= theta) .* dist(m,1);

y = cconv(a0, x0, m);

lambda = 0.1;

%% Calculate phi over the simplex
gsamps = linspace(-0.25,1.25, 1e2);        % coordinates on grid space

% Transformation from shift space to grid space -- x = Cu + d
C = [-.5 .5; -sqrt(3)/2 -sqrt(3)/2];
Cinv = inv(C);
d = [0.5; 0.5 + sqrt(3)/4];

% Shifts
s = 4;  s = -ceil(p/s):ceil(p/s);
s = s(randperm(numel(s), 3));
s = [0 s(s~=0)];  s = s(1:3);

u = [a0; zeros(m-p,1)];
v = circshift(u, s(3));  v = v(1:p);  v = v/norm(v) - a0;
u = circshift(u, s(2));  u = u(1:p);  u = u/norm(u) - a0;

phi_g = repmat({NaN(numel(gsamps),1)}, [1 numel(gsamps)]);
parfor i = 1:numel(gsamps)
    phi = phi_fista(huber(lambda));
    for j = 1:numel(gsamps)
        uv = Cinv*([gsamps(i); gsamps(j)] - d); %#ok<MINV>
        
        a = a0 + uv(1)*u + uv(2)*v; 
        a = a/norm(a);
        [phi, phi_g{i}(j)] = evaluate(phi, y, a);
    end
end
phi_g = cell2mat(phi_g);
phidelta = 8 - min(phi_g(:));
phi_g = phi_g + phidelta;

% Plot phi surface
plotcontour = true; %#ok<*NASGU>
clf; hold off;
surf_simplex; 

%% Trajectories over hemisphere
maxit = 500;

% Fix initial point
uv = [1; 1]/3;
a = a0 + uv(1)*u + uv(2)*v; 
a = a/norm(a);

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

xypath = repmat({[C*uv + d NaN(2,maxit)]}, [numel(solvers) 1]);
costs = repmat({[phi{1}.costs(end) NaN(1, maxit)]}, [numel(solvers) 1]);
%debug = repmat({NaN(2, maxit)}, [numel(solvers) 1]);
parfor j = 1:numel(solvers)
    for i = 2:maxit+1
        solvers{j} = iterate(solvers{j});
        uv = -[u v -solvers{j}.A]\a0;  uv = uv(1:2);
        
        xypath{j}(:,i) = C*uv + d;
        
        a = [u v]*uv + a0;  a = a/norm(a);
        [phi{j}, costs{j}(i)] = evaluate(phi{j}, y, a);
        
        %debug{j}(:,i-1) = [norm(a) dot(a,solvers{j}.A)];
    end
end

%% Plotting
plotcontour = true; clf; hold off;

lgd = {'PALM', 'iPALM', 'ADMM'};
colors = [0 0.5 0; 1 0 1; 1 .5 .3];
sym = {'x', 'o', '^'};
pit = ceil(maxit*0.7);
pidxs = 1:20:pit+1;

if plotcontour
for j = 1:numel(solvers)
    colormapline(...
        xypath{j}(1,1:pit), xypath{j}(2,1:pit), costs{j}(1:pit)+phidelta+1, ...
        min(log10(linspace(1,10,pit)),1)'*colors(j,:));  hold on; 
    
    plot3(...
        xypath{j}(1,pidxs), xypath{j}(2,pidxs), costs{j}(pidxs)+phidelta+1, ...
        sym{j}, 'LineWidth', 1, 'Color', colors(j,:)); %#ok<*UNRCH>
end
end

h = cell(1,2);
for j = 1:numel(solvers)
    h{1} = [h{1}, colormapline(...
        xypath{j}(1,1:pit), xypath{j}(2,1:pit), [], ...
        min(log10(linspace(1,10,pit)),1)'*colors(j,:) ...
    )];  hold on; 
    
    h{2} = [h{2}, plot(xypath{j}(1,pidxs), xypath{j}(2,pidxs), ...
        sym{j}, 'LineWidth', 1, 'Color', colors(j,:) ...
    )];
end

surf_simplex;

legend(h{2}, lgd(1:numel(solvers)));

xlim([gsamps(1) gsamps(end)]);  
ylim([gsamps(1) gsamps(end)]);

