function [A, X, A_, X_, wgts] = center(A, X, A_, X_)
if nargin < 3;          A_ = [];                end
if nargin < 4;          X_ = [];                end

for k = 1:numel(A)
    p = size(A{k});  p = p(1:2);
    p0 = floor((p+2)/3);

    wgts = cconvfft2(ones(p0), sum(abs(A{k}).^2,3));
    wgts = circshift(wgts,-floor(p0/2));

    [~,tmp] = max(wgts(:));
    [i,j] = ind2sub(p, tmp);
    i = i-ceil((p(1)+1)/2);  j = j-ceil((p(2)+1)/2);

    tmp = zeros(2*p);
    tmp(1:p(1), 1:p(2)) = A{k};
    tmp = circshift(tmp, -[i j]);
    A{k} = tmp(1:p(1), 1:p(2));
    A{k} = A{k}/norm(A{k},'fro');

    if ~isempty(A_)
        tmp(1:p(1), 1:p(2)) = A_{k};
        tmp = circshift(tmp, -[i j]);
        A_{k} = tmp(1:p(1), 1:p(2));
        A_{k} = A_{k}/norm(A_{k},'fro');
    end
    for n = 1:size(X,2)
        X{k,l} = circshift(X{k,l},-[i j]);

        if ~isempty(X_);  X_{k,l} = circshift(X_{k,l},[i j]);  end
    end
end
end

