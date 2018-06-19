classdef huber < handle & matlab.mixin.SetGet & matlab.mixin.Copyable %#ok<*PROPLC>
%HUBER  Huber function

properties
    weights = 1;      % e.g. lambda, or array for reweighting.
    mu = 1e-12;
    xpos = false;
end

methods
function obj = huber(weights, xpos, mu)
%HUBER  Constructor
    if nargin >= 1 && ~isempty(weights)
        obj.weights = weights;
    end
    if nargin >= 2 && ~isempty(xpos)
        obj.xpos = xpos;
    end
    if nargin >= 3 && ~isempty(mu)
        obj.mu = mu;
    end
end

function hx = value(obj, x)
%VALUE_  Get value for Huber function
    mu = obj.mu;

    if obj.xpos && min(x(:)) < 0
        hx = Inf;
    else
        leq = abs(x(:)) <= mu;
        tmp = NaN(numel(x),1);
        tmp(leq) = x(leq).^2/2/mu;
        tmp(~leq) = abs(x(~leq)) - mu/2;
        hx = sum(obj.weights(:) .* tmp);
    end
end

function proxh = prox(obj, x, t)
%PROX_  Compute proximal mapping for Huber function
    w_t = obj.weights/t;
    if isscalar(w_t);  w_t = w_t * ones(size(x));  end

    leq = abs(x) <= w_t + obj.mu;
    proxh = zeros(size(x));
    proxh(leq) = x(leq)./(1+w_t(leq)/obj.mu);
    proxh(~leq) = x(~leq) - w_t(~leq).*sign(x(~leq));

    if obj.xpos;  proxh = max(proxh,0);  end
end
end
end

%{
function obj = compile_params(obj)
%COMPILE_PARAMS  Recreates the function handles to update parameters
    assert(min(obj.weights(:)) >= 0, 'Weights must be nonnegative.');

    obj.value = @(x) value_(obj, x);
    obj.prox = @(x, t) prox_(obj, x, t);
    obj.diffsubg = @(x, y) diffsubg_(obj, x, y);
end

function [ eps ] = diffsubg_(obj, x, y)
%DIFFSUBG_  Difference to the subgradient of the Huber function at x.
%  In this case the Huber loss function is smooth so only need to subtract
%  from the gradient.
    w = obj.weights;
    leq = abs(x) <= obj.mu;
    subg = zeros(size(x));

    subg(leq) = x(leq)/obj.mu;
    subg(~leq) = sign(x(~leq));
    if obj.xpos;  subg(x<0) = 0;  end

    if size(w) ~= size(x);  w = reshape(w, size(x));  end
    eps = y - w.*subg;
end
%}