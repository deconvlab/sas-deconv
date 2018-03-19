classdef admm_sbd
properties
    A; X; b;                % variables
    A_; X_; Z_;             % slack variables
    WA_; WX; WZ_;           % multipliers
    
    rhoA = 1e2;  
    rhoX = 1e2;  
    rhoZ = 1e1;
    
    reg;                    % regularizer
    
    cost;
    it;
end

properties (Access = private)
    Y_;
end
    
methods
function o = admm_sbd(Y, p, reg)
    m = size(Y);
    o.Y_ = fft2(Y);
    
    tmp = [randi(m(1)) randi(m(2))];
    tmp = {mod(tmp(1)+(1:p(1)), m(1))+1 mod(tmp(2)+(1:p(2)), m(1))+1};
    
    o.A = Y(tmp{1}, tmp{2});
    %o.A = randn(p);
    o.A = o.A/norm(o.A,'fro');
    o.A_ = fft2(o.A, m(1), m(2));
    o.WA_ = zeros(m);
    
    o.X = ones(m);
    o.X_ = fft2(o.X);
    o.WX = zeros(m);
    
    o.b = 0;
    
    o.Z_ = zeros(m);
    o.WZ_ = zeros(m);
    
    o.reg = reg;
    o.it = 0;
end

function o = iterate(o)
    m = size(o.Y_);  p = size(o.A);
    
    % X_ subproblem
    tmp = o.Z_ + o.Y_;
    tmp(1) = tmp(1) - numel(o.Y_) * o.b;
    tmp = o.rhoZ * tmp - o.WZ_;
    o.X_ = fft2(o.rhoX * o.X - o.WX) + conj(o.A_) .* tmp;
    o.X_ = o.X_./(o.rhoX + o.rhoZ * abs(o.A_).^2);
    
    % X subproblem
    ifftX_ = real(ifft2(o.X_));
    %o.X = o.rhoX * o.reg.prox(ifftX_ + o.WX/o.rhoX, o.rhoX);
    o.X = ifftX_ + o.WX/o.rhoX;
    o.X = sign(o.X) .* max(abs(o.X) - o.reg.weights/o.rhoX, 0);
    
    % A_ subproblem
    o.A_ = o.rhoA * fft2(o.A, m(1), m(2)) - o.WA_ + conj(o.X_) .* tmp;
    o.A_ = o.A_./(o.rhoA + o.rhoZ * abs(o.X_).^2);
    
    % A subproblem
    tmp = real(ifft2(o.A_ + o.WA_/o.rhoA));
    tmp = tmp(1:p(1), 1:p(2));
    o.A = tmp/norm(tmp,'fro');
    
    % b subproblem
    tmp = o.Y_ + o.Z_ - o.A_ .* o.X_;
    o.b = sum(sum(real(ifft2(tmp - o.WZ_/o.rhoZ))))/numel(o.Y_);
    
    % Z_ subproblem
    tmp = o.A_ .* o.X_ - o.Y_;
    tmp(1) = tmp(1) + numel(o.Y_) * o.b;
    o.Z_ = (o.rhoZ * tmp + o.WZ_)/(1+o.rhoZ);
    
    % Update multipliers
    o.WA_ = o.WA_ + o.rhoA * (o.A_ - fft2(o.A, m(1), m(2)));
    o.WX = o.WX + o.rhoX * (ifftX_ - o.X);
    o.WZ_ = o.WZ_ + o.rhoZ * (tmp - o.Z_);
    
    % Update rho
    o.it = o.it + 1;
    if mod(o.it, 500) == 0
        o.rhoA = min(o.rhoA*2, 1e4);
        o.rhoX = min(o.rhoX*1.2, 1e4);
    end
    
    % Costs
    o.cost = zeros(4,1);
    o.cost(1) = norm(ifft2(o.Z_), 'fro')^2/2 + o.reg.weights * norm(o.X(:),1);
    %o.reg.value(o.X);
    tmp = o.A_ .* o.X_ - o.Y_ - o.Z_;
    tmp(1) = tmp(1) + numel(o.Y_) * o.b;
    o.cost(2) = norm(tmp, 'fro');
    o.cost(3) = norm(o.A_ - fft2(o.A, m(1), m(2)), 'fro');
    o.cost(4) = norm(ifftX_ - o.X, 'fro');
end
end
end