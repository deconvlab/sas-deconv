function [solver] = reweight(solver, lambda, eps, weighting)
%REWEIGHT  Apply a reweighting rule to the iPALM iterator.

if nargin < 4 || isempty(weighting) || strcmp(weighting, 'log')
    weighting = @(X) lambda./(abs(X) + eps);
elseif strcmp(weighting, 'atan')
    weighting = @(X) lambda./(abs(X).^2 + eps^2);
end

solver = reset(solver, solver.A, solver.X, solver.b);

for i = 1:numel(solver.f)
    f = solver.f{i};
    f.weights = weighting(solver.X{i}); 
    f = compile_params(f);
    solver.f{i} = f;
end
end

