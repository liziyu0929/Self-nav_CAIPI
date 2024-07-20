% Script for structured low-rank reconstruction of 2D phase map from Self-nav CAIPI
%
% Ziyu Li, Wenchuan Wu, University of Oxford, 2024


clc, clear;

%% data and parameters

addpath(genpath('utils'));

load dwi_1.22.mat
load gre_calib_1.22.mat
smoothing = 1; % apply Gaussian smoothing for output phase maps

[nx, ny, nz, nch] = size(imrecon);
nshot = 12;
im_ref = sos(imrecon);

k_data_allshot = zeros(nx, ny, nz, nch, nshot);
mask_allshot = zeros(nx, ny, nz, nshot);

load samp_180.mat
mask2show = sum(mask_allshot, 4);
figure, imshow(squeeze(mask2show(1,:,:))', [])

load ref_phs_map_1.22.mat
phs_map_ref = phs_map;

%% simulate phase-corrupted multi-shot data

for ishot = 1:nshot
    data_full = fft3c(bsxfun(@times, imrecon, phs_map_ref(:,:,ishot)));
    mask_tmp = squeeze(mask_allshot(:,:,:,ishot));
    k_data_allshot(:,:,:,:,ishot) = bsxfun(@times, data_full, mask_tmp);
end

% reference fully sampled data
kz = 7;
data_ref = fft3c(imrecon);
imref = sos(ifft2c(data_ref(:,:,kz,:)));
% figure, imshow(imref, [0 max(imref(:))*0.8])

%% SPIRiT-SLR for central kz plane
kSize = [5,5]; CalibTyk = 0.01; ReconTyk = 2e-5;
kernel = calibSPIRiT(squeeze(data_calib(:,:,7,:)), kSize, nch, CalibTyk);
GOP = SPIRiT(kernel, 'fft', [nx, ny]);

shot_full = nshot; % last kz=0 traversal shot
y_shot_full = squeeze(k_data_allshot(:,:,kz,:,shot_full));
[res_cg_x, RESVEC] = cgSPIRiT(y_shot_full, GOP, 30, ReconTyk, y_shot_full);
im_shot_full = ifft2c(res_cg_x);

y = squeeze(k_data_allshot(:,:,kz,:,:));
niter = 50; 
samp = repmat(reshape(mask_allshot(:,:,kz,:), nx, ny, 1, nshot), [1,1,nch,1]);

im_kz0_ms = spirit_slr_admm(y, niter, samp, GOP, y, abs(im_shot_full));

% combine multi-coil data
smap_central_kz = calc_sensemap([nx,ny,1,nch],data_calib(:,:,7,:));
smap = squeeze(smap_central_kz);
im_comb = sum(bsxfun(@times, im_kz0_ms, conj(smap)),3);

% Gaussian smoothing
if smoothing
    wdk2 = hamming(32);
    wdk2 = wdk2*wdk2';
    wdk2 = zpad(wdk2,[nx,ny]);

    kdata_lowrank_tmp = bsxfun(@times, fft2c(im_kz0_ms), wdk2);
    im_comb = sum(bsxfun(@times, ifft2c(kdata_lowrank_tmp), conj(smap)),3);
end

% extract phase
phs_est=squeeze(im_comb./abs(im_comb));
phs_est(find(isnan(phs_est))) = 1;

%% save results
phs_map = phs_est;
phs_est_diff = phs_est ./ squeeze(phs_est(:,:,nshot));
figure,montage(permute(angle(phs_est_diff(:,end:-1:1,:)),[2,1,3]),'Size',[3,4],'DisplayRange',[])
im = im_comb;
save('slr_phs_map_1.22.mat','phs_map','im')
