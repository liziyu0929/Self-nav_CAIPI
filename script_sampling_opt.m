% Script for sampling optimization for Self-nav CAIPI
%
% Ziyu Li, Wenchuan Wu, University of Oxford, 2024

clc, clear;

addpath(genpath('utils'));

%% specify parameters
nx = 180; ny = 180; nz = 12; nch = 8; nshot = 12;

k_data_allshot = zeros(nx, ny, nz, nch, nshot);
mask_allshot = zeros(nx, ny, nz, nshot);

kz_cover = 6; % kz_coverage
p = [0:1:kz_cover, kz_cover-1:-1:1];  % one period of sampling

N = 1e8; 
opt_para = zeros(3, nshot);
tmp_para = zeros(3, nshot);
metrics = 1e6; 
kcntr_dist = 0; ovlp = 0; uni = 0;

% one shot covering the central kz=0
mask_allshot(:, 2:3:end, 7, 12) = 1;

kcntr_dist_allshot = 0;
for ishot = 1:11
    disp(['searching shot' num2str(ishot)]);
    metrics = 1e6; 
    if ishot == 1 % no need to calculate k-space gaps for the first CAIPI shot
        for kysft = 1:3
            for psft = 1:11
                kzsft = 1;
                flag = 1; kcntr_dist_tmp = 1e6;
                shot_tmp = samp_shot(nx,ny,nz,kysft,kzsft,p,psft,kz_cover);
                itsct = find(shot_tmp(1,:,7));
                for ii = 1:length(itsct)
                    if (length(itsct) > 6) && (abs(itsct(ii) - (ny/2+1)) < 15) || (length(itsct) <= 6) && (abs(itsct(ii) - (ny/2+1)) < 3)
                        kcntr_dist_tmp = min(kcntr_dist_tmp, abs(itsct(ii) - (ny/2+1)));
                        flag = 0;
                    end
                    if ii > 1
                        a(ii-1)=itsct(ii)-itsct(ii-1);
                    end
                end
                if flag
                    continue
                end
                mask_allshot_tmp = mask_allshot;
                mask_allshot_tmp(:,:,:,ishot) = shot_tmp;
                
                maskall = sum(mask_allshot_tmp, 4);
                mask_yz = squeeze(maskall(1,:,:));
                
                idx = find(mask_yz>1);
                ovlp_tmp = sum(mask_yz(idx));
                metrics_tmp = ovlp_tmp + kcntr_dist_tmp;
                if metrics_tmp < metrics
                    metrics = metrics_tmp; 
                    kcntr_dist = kcntr_dist_tmp + kcntr_dist_allshot;
                    ovlp = ovlp_tmp;
                    disp(['current metrics:' num2str(metrics) ' ovlp:' ...
                        num2str(ovlp) ' kcntr_dist:' num2str(kcntr_dist)])
                    best = shot_tmp;
                    opt_para(:,ishot) = [kysft,kzsft,psft];
                    best_kcntr_dist = kcntr_dist_tmp;
                end
            end
        end
        mask_allshot(:, :,:, 1) = best;
        kcntr_dist_allshot = kcntr_dist_allshot + best_kcntr_dist;
    
    else
    for kysft = 1:3
        for kzsft = 1:6
            for psft = 1:11
                flag = 1; kcntr_dist_tmp = 1e6;
                shot_tmp = samp_shot(nx,ny,nz,kysft,kzsft,p,psft,kz_cover);
                itsct = find(shot_tmp(1,:,7));
                for ii = 1:length(itsct)
                    if (length(itsct) > 6) && (abs(itsct(ii) - (ny/2+1)) < 15) || (length(itsct) <= 6) && (abs(itsct(ii) - (ny/2+1)) < 3)
                        kcntr_dist_tmp = min(kcntr_dist_tmp, abs(itsct(ii) - (ny/2+1)));
                        flag = 0;
                    end
                    if ii > 1
                        a(ii-1)=itsct(ii)-itsct(ii-1);
                    end
                end
                if flag
                    continue
                end
                mask_allshot_tmp = mask_allshot;
                mask_allshot_tmp(:,:,:,ishot) = shot_tmp;
                
                maskall = sum(mask_allshot_tmp, 4);
                mask_yz = squeeze(maskall(1,:,:));
                if ishot > 1 
                    ker = ones(3,3); 
                    conv = conv2(mask_yz, ker, 'same'); 
                    gap3 = length(find(conv==0));
                    gap_tmp = gap3;
                else
                    gap_tmp = 0;
                end
                
                idx = find(mask_yz>1);
                ovlp_tmp = sum(mask_yz(idx));
                metrics_tmp = ovlp_tmp + kcntr_dist_tmp + gap_tmp + kcntr_dist_allshot;
                if metrics_tmp < metrics
                    metrics = metrics_tmp; 
                    kcntr_dist = kcntr_dist_tmp + kcntr_dist_allshot;
                    ovlp = ovlp_tmp; uni = gap_tmp;
                    disp(['current metrics:' num2str(metrics) ' gap:' num2str(uni) ' ovlp:' ...
                        num2str(ovlp) ' kcntr_dist:' num2str(kcntr_dist)])
                    best = shot_tmp;
                    opt_para(:,ishot) = [kysft,kzsft,psft];
                    best_kcntr_dist = kcntr_dist_tmp;
                end
            end
        end
    end
    mask_allshot(:,:,:,ishot) = best;
    kcntr_dist_allshot = kcntr_dist_allshot + best_kcntr_dist;
    end
end

mask = sum(mask_allshot, 4);
figure, imshow(squeeze(mask(1,:,:))', [])
