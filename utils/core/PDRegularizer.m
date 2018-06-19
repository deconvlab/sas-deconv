classdef PDRegularizer < handle & matlab.mixin.SetGet & matlab.mixin.Copyable %#ok<*PROPLC>
%PDREGULARIZER  A primal-dual regularizer of the form r(X) = g(DX)

properties
    g;  weights;            % function with a prox; weights shared with g
    D;                      % either a matrix or function handles for D and D*

    x;  y;                  % primal / dual variable
    xintv = [-Inf Inf];
    stepsz = 1e-1;          % stepsize on primal variable
    
    niter = 10;             % # of iterations to take when calculating prox
    cost;  it = 0;
end

methods
function o = PDRegularizer(g, D, weights, stepsz, xintv)
    o.g = g;

    if nargin < 2 || isempty(D)
        o.D = struct('funh', @(x) x, 'adjh', @(y) y, 'stepsz', 1);
    elseif isnumeric(D) && ismatrix(D)
        o.D = struct('mtx', D, 'funh', @(x) D*x, 'adjh', @(y) D'*y, 'stepsz', 1/svds(D,1)^2);
    else
        o.D = D;
    end

    if nargin < 3 || isempty(weights)
        o.weights = g.weights;
    else
        o.weights = weights;
    end
    
    if nargin >= 4 && ~isempty(stepsz);  o.stepsz = stepsz; end
    if nargin >= 5 && ~isempty(xintv);   o.xintv = xintv;   end
end

function o = reset(o, x, y)
    if nargin >= 2 && ~isempty(x);      o.x = x;  end
    if nargin >= 3 && ~isempty(y);      o.y = y;  end
    if isempty(o.y);  o.y = o.D.funh(x);  end

    o.it = 0;       o.cost = [];
end

function hx = value(o, x)
    o.g.weights = o.weights;
    hx = o.g.value(o.D.funh(x));
end

function o = step(o, x_, t)
%STEP  Take a step using primal-dual
    o.g.weights = o.weights;
    stepsz = 0.99 * [o.stepsz, o.D.stepsz/o.stepsz];

    % update x
    r = o.x - stepsz(1)/t * o.D.adjh(o.y);
    xnew = (r + 2*stepsz(1)*x_)/(1+2*stepsz(1));
    o.x = 2 * xnew - o.x;
    o.x = max(o.x, o.xintv(1));
    o.x = min(o.x, o.xintv(2));
    
    % update y
    tmp = stepsz(2)/t;
    r = o.y + tmp * o.D.funh(o.x);
    o.y = r - tmp * o.g.prox(r/tmp, tmp);

    o.cost = [norm(o.x(:)-x_(:))^2/2 o.g.value(o.D.funh(o.x))/t];
    o.it = o.it + 1;
end

function [out, costs] = prox(o, x_, t, niter)
    if nargin < 3 || isempty(t);  t = 1;  end
    if isempty(o.x);  o.reset(x_);  end

    if nargin < 4 || isempty(niter);  niter = o.niter;  end
    costs = NaN(2, niter);
    for it = 1:niter
        o = o.step(x_, t);
        costs(:,it) = o.cost;
    end
    out = o.x;
    costs = costs';
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

function o = setparams(o, theta, alph, stepsz)
%SETPARAMS  A way for users to set the APD parameters. Usually there is no
%   convergence issues even with default settings so rarely needs to be used.
%
    if nargin >= 2 && ~isempty(theta);  o.theta = theta;  end
    if nargin >= 3 && ~isempty(alph);   o.alph = alph;  end
    if nargin >= 4 && ~isempty(stepsz); o.stepsz = stepsz;  end
end

%}