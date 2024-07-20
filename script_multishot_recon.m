
clc, clear

%% data and parameters

addpath(genpath('utils'));

load dwi_1.22.mat
load gre_calib_1.22.mat

im_ref = sos(imrecon);

[nx, ny, nz, nch] = size(imrecon); nshot = 12;

load ref_phs_map_1.22.mat
phs_map_ref = phs_map;

k_data_allshot = zeros(nx, ny, nz, nch, nshot);
mask_allshot = zeros(nx, ny, nz, nshot);

% fully sampled reference data
figure,...
    montage(permute(im_ref(:,end:-1:1,:),[2,1,3]),...
    'Size',[3,4],'DisplayRange',[0 max(im_ref(:))*0.4])


load samp_180.mat

mask2show = sum(mask_allshot, 4);
figure, imshow(squeeze(mask2show(1,:,:))', [])

%% simulate phase-corrupted multi-shot data

for ishot = 1:nshot
    data_full = fft3c(bsxfun(@times, imrecon, phs_map_ref(:,:,ishot)));
    mask_tmp = squeeze(mask_allshot(:,:,:,ishot));
    k_data_allshot(:,:,:,:,ishot) = bsxfun(@times, data_full, mask_tmp);
end


%% SPIRiT phase-corrected multi-shot recon
load slr_phs_map_1.22.mat
phs_map = repmat(phs_map, [1 1 1 nch]);
phs_map = squeeze(phs_map);

% SPIRiT
kSize = [5,5]; % SPIRiT kernel size
lambda = 2e-4; % tiknov
alpha = 1; % SPIRiT regularization
nIterCG = 30;
CalibTyk = 0.1;

kdata_ds_1 = ifft1c(k_data_allshot,1);
data_calib_fro = ifft1c(data_calib,1);

nz_eff = size(kdata_ds_1, 3);
k_recon = zeros(nx,ny,nz,nch);
samp = mask2samp(mask_allshot); % from k-space mask to sampling matrix

parfor ix = 1:nx
    disp(num2str(ix))
    sz = [ny,nz_eff];
    kdata_x = squeeze(kdata_ds_1(ix,:,:,:,:));
    data_phs_x = squeeze(phs_map(ix,:,:,:));
    
    DATA_x_calib = squeeze(data_calib_fro(ix,:,:,:));

    kernel = calibSPIRiT(DATA_x_calib, kSize, nch, CalibTyk);
    GOP = SPIRiT(kernel, 'fft', [ny, nz_eff]);
    [resSPIRiT] = pcgSPIRiT_multishot_phaseCorrection(...
                    kdata_x, data_phs_x, samp, GOP, nIterCG, sz, nch, 1, lambda, alpha);

    k_recon(ix,:,:,:)=resSPIRiT;
     
end
imrecon=sos(ifft3c(fft1c(k_recon,1)));

% reconstruction results
figure,...
    montage(permute(imrecon(:,end:-1:1,:),[2,1,3]),...
    'Size',[3,4],'DisplayRange',[0 max(imrecon(:))*0.4])

% residual map
res = im_ref./max(im_ref(:)) - imrecon ./ max(im_ref(:));
figure,...
    montage(permute(res(:,end:-1:1,:),[2,1,3]),...
    'Size',[3,4],'DisplayRange',[-0.1 0.1])