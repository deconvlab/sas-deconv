function imgupdate(Y, obj, A0, X0)
A = obj.A{1};  X = obj.X{1}; it = obj.it;
p = size(A);

if nargin >= 3
    p0 = size(A0);
    tmp = zeros(p);
    tmp(floor(p(1)/2)+((-ceil(p0(1)/2)+1):floor(p0(1)/2)), ...
        floor(p(2)/2)+((-ceil(p0(2)/2)+1):floor(p0(2)/2))) = A0;
    subplot(321);  imagesc(abs(Y)); title(['Iteration ' num2str(it)]);
    subplot(323);  imagesc(abs(tmp));
    subplot(325);  imagesc(abs(circshift(X0, floor(p0/2))));

    subplot(322);  imagesc(abs(cconvfft2(A, X)));
    subplot(324);  imagesc(abs(A));
    subplot(326);  imagesc(abs(circshift(X, floor(p/2))));
else
    subplot(221);  imagesc(abs(Y)); title(['Iteration ' num2str(it)]);
    subplot(222);  imagesc(abs(cconvfft2(A, X)));
    subplot(223);  imagesc(abs(A));
    subplot(224);  imagesc(abs(circshift(X, floor(p/2))));
end
drawnow;
end

