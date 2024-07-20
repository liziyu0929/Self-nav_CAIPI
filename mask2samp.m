% from mask to sampling matrix
%
% Ziyu Li, Wenchuan Wu, University of Oxford, 2024

function samp_order = mask2samp(mask_allshot)
[nx, ny, nz, nshot] = size(mask_allshot);
% samp_order = zeros(ceil(ny/Ry),2,1);

for ii = 1 : nshot
    idx_samp = find(squeeze(mask_allshot(1,:,:,ii))==1);
    [ky, kz] = ind2sub(size(squeeze(mask_allshot(1,:,:,ii))), idx_samp);
    [ky_s, idx] = sort(ky);
    kz_s = kz(idx);
    samp_order(:, 1, ii) = ky_s; % ky sampling
    samp_order(:, 2, ii) = kz_s; % kz sampling
end
end
    