function samples = sample_hemi3(nr, nt)
%SAMPLE_HEMI3  Sample points across the positive 3-hemisphere.
[t, r] = meshgrid(linspace(0,2*pi,nt),linspace(0,1,nr));
samples = zeros([nr nt 3]);

samples(:,:,1) = r.*cos(t);  
samples(:,:,2) = r.*sin(t);
samples(:,:,3) = sqrt(max(1-samples(:,:,1).^2-samples(:,:,2).^2, 0));
end

function samples = sample_hemi3_rings(samp_rings, samps_bdy) %#ok<*DEFNU>
%SAMPLE_HEMI3  Sample points across the positive 3-hemisphere.
%  
%   samples = sample_hemi3(samp_layers, samps_bdy)
%   
%   Arguments:
%     samp_rings:  number of rings to sample between the origin and the
%       boundary (x^2 + y^2 = 1). The total number of rings is samp_layers 
%       + 2 and the radius of each ring is uniform between 0 and 1.
%
%     samps_bdy:  number of samples to take on the boundary. The number of
%     samples for each ring is (roughly) proportional to its radius.
%
%
%   Output:
%     samples:  a cell array of length (samp_rings + 2) containing samples
%       at each point. Each cell element contains an array of size
%       (n_samps x 3).
%
%

samples = cell(samp_rings+2,1);
samples{1} = [0 0 1];

r = linspace(0,1,samp_rings+2);
n = max(ceil(r*samps_bdy), 3);

for i = 2:samp_rings+2
    theta = linspace(0, 2*pi, n(i)+1)';
    theta = theta(1:end-1);
    samples{i} = [r(i)*cos(theta) r(i)*sin(theta) sqrt(ones(n(i),1)-r(i)^2)];
end
end

