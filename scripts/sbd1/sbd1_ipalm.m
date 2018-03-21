classdef sbd1_ipalm
properties
    A
    X
    iterator
end
methods
function o = sbd1_ipalm(y, p, lambda, xpos, getbias, a0)
    %SBD1_IPALM  SBD in 1D using iPALM

    % Cost function + some defs
    m = numel(y);  s = sphereops;  y = y(:);
    H.value = @(A, X, b, c) Hval(A{:}, X{:}, b{:}, y, c);


    % Set up A
    if nargin < 6;  a0 = y(mod(randi(m)+(1:p(1)), m(1))+1);  end
    a0 = s.proj(a0(:));
    H.gradA = {@(A, X, b, c) grada(A{:}, X{:}, b{:}, y, c)};
    tA = {@stepsza};


    % Set up X
    H.gradX = {@(A, X, b, c) gradx(A{:}, X{:}, b{:}, y, c)};
    tX = {@stepszx};
    f = {huber(lambda, xpos)};


    % Set up b
    H.gradb = {@(A, X, b, c) gradb(A{:}, X{:}, b{:}, y, getbias, c)};
    tb = {@(~,~,~,c) stepszb([],[],[], y, getbias, c)};
    if getbias;  b0 = median(y);  else;  b0 = 0;  end

    % Set up the object properties
    o.iterator = ipalm(H, f, {a0}, {rand(m,1)}, {b0}, tA, tX, tb);
    o.A = o.iterator.A{1};
    o.X = o.iterator.X{1};
end

function o = iterate(o)
    o.iterator = iterate(o.iterator);
    o.A = o.iterator.A{1};
    o.X = o.iterator.X{1};
end
end
end

function [v, cache] = Hval(a, x, b, y, cache)
if nargin < 5 || isempty(cache)
    v = norm(cconv(a, x, numel(x)) + b - y)^2/2;
else
    v = cache.Hcost;
end
end

function [g, cache] = grada(a, x, b, y, cache)
xhat = fft(x);  
g = ifft(conj(xhat) .* (fft(a,numel(y)) .* xhat + fft(b-y)));
g = real(g(1:numel(a)));

if nargin >= 5 && ~isempty(cache)
    cache.tA = 1/max(abs(xhat(:)))^2;
end
end

function [t, cache] = stepsza(~, ~, ~, cache)
t = cache.tA;
end

function [g, cache] = gradx(a, x, b, y, cache)
ahat = fft(a, numel(y));  
g = real(ifft(conj(ahat) .* (ahat .* fft(x) + fft(b-y))));

if nargin >= 5 && ~isempty(cache)
    cache.tX = 1/max(abs(ahat(:)))^2;
end
end

function [t, cache] = stepszx(~, ~, ~, cache)
t = cache.tX;
end

function [g, cache] = gradb(a, x, b, y, getbias, cache)
r = cconv(a, x, numel(x)) + b - y;  
if getbias;  g = sum(r(:));  else;  g = 0;  end

if nargin >= 6 && ~isempty(cache)
    cache.Hcost = norm(r)^2/2;
end
end

function [t, cache] = stepszb(~, ~, ~, y, getbias, cache)
if getbias
    t = 1/(2 * numel(y));  
else
    t = 0;  
end
end
