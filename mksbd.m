function iterator = mksbd(Y, p, lambda, xpos, getbias)
%MKSBD  Make an iterator for solving SBD
%  Creates an iPALM iterator for solving SBD by supplying a smooth
%  square-error term H, the nonsmooth Huber terms and the corresponding
%  descent (gradient / prox quantities).
%
%  iterator = mkcdl(Y, p, lambda, xpos, getbias)
%    returns the iterator. Provide as input arguments:
%       Y:  [ARRAY DOUBLE].  A cell array containing the observation.
%       p:  [2 INT].  The size of the recovered kernel.
%       lambda:  [DOUBLE >0].  The sparsity tradeoff parameter.
%       xpos:  bool.  Set to TRUE to ensure X is nonnegative.
%       getbias: bool.  Set to TRUE to estimate a constant bias for Y.
%

% Cost function + some defs
m = size(Y);  ob = obops;  
H.value = @(A, X, b, c) Hval(A{:}, X{:}, b{:}, Y, c);


% Set up A
tmp = [randi(m(1)) randi(m(2))];
tmp = {mod(tmp(1)+(1:p(1)), m(1))+1 mod(tmp(2)+(1:p(2)), m(1))+1};
A0 = ob.proj(Y(tmp{1}, tmp{2}));

H.gradA = {@(A, X, b, c) gradA(A{:}, X{:}, b{:}, Y, c)};
tA = {@stepszA};


% Set up X
H.gradX = {@(A, X, b, c) gradX(A{:}, X{:}, b{:}, Y, c)};
tX = {@stepszX};
f = {huber(lambda, xpos)};


% Set up b
H.gradb = {@(A, X, b, c) gradb(A{:}, X{:}, b{:}, Y, getbias, c)};
tb = {@(~,~,~,c) stepszb([],[],[], Y, getbias, c)};


iterator = ipalm(H, f, {A0}, {ones(m)}, {median(Y(:))}, tA, tX, tb);
end

function [v, cache] = Hval(A, X, b, Y, cache)
if nargin < 5 || isempty(cache)
    v = norm(cconvfft2(A, X) + b - Y, 'fro')^2/2;
else
    v = cache.Hcost;
end
end

function [g, cache] = gradA(A, X, b, Y, cache)
Xhat = fft2(X);  
g = ifft2(conj(Xhat) .* (fft2(A,size(Y,1),size(Y,2)) .* Xhat + fft2(b-Y)));
g = real(g(1:size(A,1), 1:size(A,2)));

if nargin >= 5 && ~isempty(cache)
    cache.tA = 1/max(abs(Xhat(:)))^2;
end
end

function [t, cache] = stepszA(~, ~, ~, cache)
t = cache.tA;
end

function [g, cache] = gradX(A, X, b, Y, cache)
Ahat = fft2(A, size(Y,1), size(Y,2));  
g = real(ifft2(conj(Ahat) .* (Ahat .* fft2(X) + fft2(b-Y))));

if nargin >= 5 && ~isempty(cache)
    cache.tX = 1/max(abs(Ahat(:)))^2;
end
end

function [t, cache] = stepszX(~, ~, ~, cache)
t = cache.tX;
end

function [g, cache] = gradb(A, X, b, Y, getbias, cache)
R = cconvfft2(A, X) + b - Y;  
if getbias;  g = sum(R(:));  else;  g = 0;  end

if nargin >= 6 && ~isempty(cache)
    cache.Hcost = norm(R, 'fro')^2/2;
end
end

function [t, cache] = stepszb(~, ~, ~, Y, getbias, cache)
if getbias
    t = 1/(2 * numel(Y));  
else
    t = 0;  
end
end
