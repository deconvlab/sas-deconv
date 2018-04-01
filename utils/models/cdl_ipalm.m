classdef cdl_ipalm < ipalm
properties
    synthesize = @synthfun;
end  
    
methods
function o = cdl_ipalm(Y, p, K, lambda, xpos, getbias)
%CDL_IPALM Iterator for solving CDL
%  Creates an iPALM iterator for solving CDL by supplying a smooth
%  square-error term H, the nonsmooth Huber terms and the corresponding
%  descent (gradient / prox quantities).
%
%  [iterator, synthesize] = mkcdl(Y, p, K, lambda, xpos, getbias)
%    returns the iterator and a function for synthesizing a reconstruction
%    of the data using cell arrays A and X in the form
%           A:[1xK CELL],  X:[KxN CELL]
%    where K is the number of kernels and N is the number of samples in Y.
%
%    Provide as input arguments:
%       Y:  [1xN CELL].  A cell array containing the observation samples.
%       p:  [2 INT].  The size of the recovered kernel.
%       K:  [INT].  The number of kernels to recover.
%       lambda:  [DOUBLE >0].  The sparsity tradeoff parameter.
%       xpos:  bool.  Set to TRUE to ensure X is nonnegative.
%       getbias: bool.  Set to TRUE to estimate a constant bias in each Y.
%

    % Cost function + some defs
    N = numel(Y);
    ob = obops;  
    afun = @(funh, a) arrayfun(funh, a{:}, 'UniformOutput', false);
    H.value = @(A, X, b, c) Hval(A, X, b, Y, c);


    % Set up A
    A0 = cell(1, K);
    for i = 1:K
        idx = randi(N);
        m = size(Y{idx});
        tmp = [randi(m(1)) randi(m(2))];
        tmp = {mod(tmp(1)+(1:p(1)), m(1))+1 mod(tmp(2)+(1:p(2)), m(1))+1};
        A0{i} = ob.proj(Y{idx}(tmp{1}, tmp{2}));
    end

    H.gradA = afun(@(i) @(A, X, b, c) gradA(A, X, b, Y, i, c), {1:K});
    tA = afun(@(i) @(A, X, b, c) stepszA(A, X, b, i, c), {1:K});


    % Set up X
    idx = {repmat((1:K)', [1 N])  repmat((1:N), [K 1])};
    H.gradX = afun(@(k, n) @(A, X, b, c) gradX(A, X, b, Y, [k n], c), idx);
    tX = afun(@(k, n) @(A, X, b, c) stepszX(A, X, b, [k n], c), idx);
    f = afun(@(~) huber(lambda, xpos), idx(1));

    X0 = repmat({ones(m)}, [K N]);


    % Set up b
    H.gradb = afun(@(i) @(A, X, b, c) gradb(A, X, b, Y, i, getbias, c), {1:N});
    tb = afun(@(i) @(A,~,~,c) stepszb(A,[],[], Y, i, getbias, c), {1:N});

    if getbias
        b0 = afun(@(i) median(Y{i}(:)), {1:N});
    else
        b0 = afun(@(i) 0, {1:N});
    end

    o = o@ipalm(H, f, A0, X0, b0, tA, tX, tb);
end


function Yhat = reconstruct(o, n)
%RECONSTRUCT  Reconstruct observation based on CDL parameters
if nargin < 2 || isempty(n)
    Yhat = synthfun(o.A, o.X, o.b);
else
    Yhat = synthfun(o.A, o.X(:,n), o.b(n));
end
end

end
end


function Yest = synthfun(A, X, b, invfft)
[K, N] = size(X);
if nargin < 3 || isempty(b);   afun(@(i) 0, {ones(1,N)});     end
if nargin < 4 || isempty(invfft);   invfft = true;     end

Yest = cell(1, N);
for n = 1:N
    m = size(X{1, n});
    tmp = zeros(m);  tmp(1) = tmp(1) + prod(m) * b{n};
    for k = 1:K
        tmp = tmp + fft2(A{k}, m(1), m(2)) .* fft2(X{k, n});
    end
    if invfft;  tmp = real(ifft2(tmp));  end
    Yest{n} = tmp;
end

if N == 1;  Yest = Yest{1};  end
end


function [Rn, cache] = residual(A, X, b, Y, n, invfft, cache)
if nargin < 5 || isempty(invfft);   invfft = true;     end

if invfft
    Rn = Y{n};
elseif nargin == 6 && isstruct(cache)
    if ~isfield(cache, 'Yhat')
        cache = cell(size(Y));
        for k = 1:numel(Y)
            cache.Yhat{k} = fft2(Y{k});
        end
    end
    Rn = cache.Yhat{n};
else
    Rn = fft2(Y{n});
end

Rn = synthfun(A,X(:,n),b, invfft) - Rn; 
end


function [v, cache] = Hval(A, X, b, Y, cache)
if nargin < 5 || ~isfield(cache, 'Hcost')
    v = 0;
    for i = 1:numel(Y)
        R = residual(A, X, b, Y, n);
        v = v + norm(R(:))^2/2;
    end
else
    v = sum(cache.Hcost);
end
end


function [g_Ak, cache] = gradA(A, X, b, Y, k, cache)
LAk = 0; 
g_Ak = zeros(size(A{k}));
N = numel(Y);

for n = 1:N
    [Rihat, cache] = residual(A, X, b, Y, n, false, cache);
    Xinhat = fft2(X{k, n});  
    LAk = max(LAk, max(abs(Xinhat(:))));
    tmp = real(ifft2(conj(Xinhat) .* Rihat));
    g_Ak = g_Ak + tmp(1:size(A{k},1), 1:size(A{k},2));
end
    
if nargin >= 6 && isstruct(cache)
    if ~isfield(cache, 'LA');  cache.LA = zeros(size(A));  end
    cache.LA(k) = sqrt(N)*LAk^2;
end
end


function [tAk, cache] = stepszA(~, ~, ~, k, cache)
tAk = 0.99/cache.LA(k);
end


function [g, cache] = gradX(A, X, b, Y, idx, cache)
[Rhat, cache] = residual(A, X, b, Y, idx(2), false, cache);
Aihat = fft2(A{idx(1)}, size(Y{idx(2)},1), size(Y{idx(2)},2));
g = real(ifft2(conj(Aihat) .* Rhat));

if nargin >= 5 && ~isempty(cache)
    if ~isfield(cache, 'LX');  cache.LX = zeros(size(X));  end
    cache.LX(idx(1), idx(2)) = max(abs(Aihat(:)))^2;
end
end


function [t, cache] = stepszX(~, ~, ~, idx, cache)
t = 0.99/cache.LX(idx(1), idx(2));
end


function [g, cache] = gradb(A, X, b, Y, n, getbias, cache)
R = residual(A, X, b, Y, n, true, cache);
if getbias
    g = sum(R(:));
else
    g = 0;  
end

if nargin >= 6 && ~isempty(cache)
    if ~isfield(cache, 'Hcost');  cache.Hcost = zeros(numel(Y),1);  end
    cache.Hcost(n) = norm(R(:))^2/2;
end
end


function [t, cache] = stepszb(A, ~, ~, Y, n, getbias, cache)
if getbias
    t = 1/(2 * numel(Y{n}) * sqrt(numel(A)));  
else
    t = 0;  
end
end