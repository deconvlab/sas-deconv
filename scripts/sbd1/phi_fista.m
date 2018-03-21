classdef phi_fista
properties
    rho;
    
    maxit = 1e4;
    it;
    
    tol = 1e-4;
    eps;
    
    costs;
    x;
end

methods
function o = phi_fista(rho)
    o.rho = rho;
end
function [o, cost] = evaluate(o, y, a)
    m = numel(y);  y = y(:);  a = a(:);
    
    ahat = fft(a,m);
    a2hat = abs(ahat).^2;
    ayhat = conj(ahat).*fft(y);
    s = 0.99/max(a2hat);
    
    if isempty(o.x)
        o.x = randn(m,1);  
    end
    
    w = fft(o.x);  xhat_ = w;  t = 1;
    o.it = 0;  repeat = true;  
    o.costs = NaN(o.maxit,1); cost = Inf;
    while repeat
        o.x = o.rho.prox(real(ifft(w - s*(a2hat.*w-ayhat))), 1/s);
        t_ = (1+sqrt(1 + 4*t^2))/2;
        xhat = fft(o.x);
        w = xhat + (t-1)/t_ * (xhat - xhat_);

        t = t_;  xhat_ = xhat;
        o.it = o.it + 1;
        
        o.costs(o.it) = norm(real(ifft(ahat.*xhat))-y)^2/2 + o.rho.value(o.x);
        o.eps = abs(cost - o.costs(o.it));
        cost = o.costs(o.it);
        repeat = (o.it < o.maxit) && (o.eps >= o.tol);
    end
    
    o.costs = o.costs(1:o.it);
end
end
end