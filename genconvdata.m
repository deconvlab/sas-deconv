function [A0, X0, Y] = genconvdata(params_in)
%GENCONVDATA  Generate simulated convolutional data
%
%  [A0, X0, Y] = genconvdata(params_in)  generates *clean* cyclic
%  convolutional data, i.e. Y = cconvfft2(A0, X0); any noise or biases
%  should be added in after, or using A0 and X0.
%  
%
%  The data is generated using the parameters supplied by the struct
%  params_in. In particular:
%
%  A0 - the kernel is generated randomly or based on simulated data, and
%  each kernel is automatically placed on the sphere. The number of
%  kernels and its sizes can be specified.
%
%  X0 - the activation map is generated based on a Bernoulli model. The
%  size of the observation, the distribution of the spikes and the
%  Bernoulli density can be specified.
%
%  N - the number of samples generated can be specified.
%
%  Please refer to the defaults in the script for reference; fields
%  specified in  params_in  override the defaults.
%
%
%  Finally, if either the number of kernels K or samples N are greater than 
%  1, then A0, X0 and Y are returned as cell arrays (otherwise they are 
%  returned simply as numerical arrays). Cell arrays are arranged in
%  matrix-vector convention:
%  
%    A0: [1 x K],  X0:[K x N],  Y:[1 x N];
%

%% Default data parameters
% Kernel
kerneltype = 'simulated';       % Either simulated STM kernels or random
simslices = [18 25];            % For STM kernels, choose the slices
simslicewin = 22;               % No need to change this unless different 
                                %   STM file is picked - full size 34

K = 2;                          % Number of kernels
p = 32;                         % Size of the (short) kernel

% Activation / observation
m = 256;                        % Observation size
theta = 3e-3;                   % Bernoulli (sparsity) coefficient
dist = @ones;                   % Distribution of activations

N = 1;                          % Number of samples

%% Process params_in
if nargin >= 1 && isstruct(params_in)
    tmp = fieldnames(params_in);
    for i = 1:numel(tmp)
        eval(sprintf('%s=getfield(params_in,''%s'');', tmp{i}, tmp{i}));
    end
end

%% Generate data
addpath('helpers');
ob = obops;

assert(p < m, 'Kernel size should be small compared to observation.');

if strcmp(kerneltype, 'random')
    A0 = randn([p p K]);
else
    assert(K<=numel(simslices), ['Number of slices chosen for ' ...
        'simulated data -- numel(simslices) -- must be at least K.']);
    
    A0 = zeros([p p K]);  win = simslicewin;
    
    load('data/Data_N_50_nDef_1_thop_-0.200.mat', 'LDoS');
    for i = 1:K
        tmp = LDoS(131+(-win:win), 131+(-win:win), simslices(i));
        A0(:,:,i) = imresize(tmp, [p p]);
    end
    clear LDoS;
end
A0 = ob.proj(A0);

X0 = double(rand([m m K N]) <= theta) .* dist([m m K N]);

Y = zeros([m m N]);
for k = 1:K
    for i = 1:N
        Y(:,:,i) = Y(:,:,i) + cconvfft2(A0(:,:,k), X0(:,:,k,i));
    end
end

if (K > 1) || (N > 1)
    A0 = squeeze(mat2cell(A0, p, p, ones(1,K)))';
    X0 = squeeze(mat2cell(X0, m, m, ones(K,1), ones(1,N)));
    Y = squeeze(mat2cell(Y, m, m, ones(1,N)))';
end
end