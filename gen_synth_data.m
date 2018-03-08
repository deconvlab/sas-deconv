%% Data parameters
% Kernel
kerneltype = 'simulated';
slices = [18 25];
slicerad = 22;   % full size 34
K = 2;
p0 = 32;

% Activation / observation
m = 256;
theta = 3e-3;
dist = ones([m m K]);

%% Generate data
ob = obops;

if strcmp(kerneltype, 'random')
    A0 = randn([p0 p0 K]);
else
    A0 = zeros([p0 p0 K]);
    load('helpers/Data_N_50_nDef_1_thop_-0.200.mat', 'LDoS');
    for i = 1:K
        tmp = LDoS(131+(-slicerad:slicerad), 131+(-slicerad:slicerad), slices(i));
        A0(:,:,i) = imresize(tmp, [p0 p0]);
    end
    clear LDoS;
end
A0 = ob.proj(A0);

X0 = double(rand([m m K]) <= theta) .* dist;

Y = zeros([m m]);
for i = 1:K
    Y = Y + cconvfft2(A0(:,:,i), X0(:,:,i));
end

A0 = squeeze(mat2cell(A0, p0, p0, ones(K,1)))';
X0 = squeeze(mat2cell(X0, m, m, ones(K,1)))';

%%
% imagesc(abs(Y));
clearvars -EXCEPT A0 X0 Y