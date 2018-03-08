function [A, X, wgts, A_, X_] = center(A, X, p0, A_, X_)
p = size(A);  p = p(1:2);

if nargin < 3 || isempty(p0);   p0 = floor((p+2)/3);    end
if nargin < 4;                  A_ = [];                end
if nargin < 5;                  X_ = [];                end
if numel(p0) == 1;              p0 = [p0 p0];           end

wgts = cconvfft2(ones(p0), sum(abs(A).^2,3));
wgts = circshift(wgts,-floor(p0/2));

[~,tmp] = max(wgts(:));
[i,j] = ind2sub(p, tmp);
i = i-ceil((p(1)+1)/2);  j = j-ceil((p(2)+1)/2);

tmp = zeros(2*p);
for k = 1:size(A,3)
    tmp(1:p(1), 1:p(2)) = A(:,:,k);
    tmp = circshift(tmp, -[i j]);
    A(:,:,k) = tmp(1:p(1), 1:p(2));
    A(:,:,k) = A(:,:,k)/norm(A(:,:,k),'fro');
    
    if ~isempty(A_)
        tmp(1:p(1), 1:p(2)) = A_(:,:,k);
        tmp = circshift(tmp, -[i j]);
        A_(:,:,k) = tmp(1:p(1), 1:p(2));
        A_(:,:,k) = A_(:,:,k)/norm(A_(:,:,k),'fro');
    end
end
X = circshift(X,-[i j]);  
if ~isempty(X_);  X_ = circshift(X_,[i j]);  end
end

