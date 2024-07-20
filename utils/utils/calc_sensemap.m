function smap = calc_sensemap(size_map, data_calib, nset, flag)

% calculate sensitivity map using ESPIRiT package
% size: size of sensitivity maps
% data_calib: calibration data
% nset: number of sets of sensitivity maps
% flag: 3D (1) or 2D (0) calculation

if nargin < 4
    flag = 0;
    nset = 1;
end

disp('calculating sensitivity maps...')
tic
smap = zeros([size_map nset]);
kSize = [5, 5];
eigThresh_1 = 0.04;
eigThresh_2 = 0.98;
data_calib_fro = ifft1c(data_calib, 3);
if flag == 0
    for isli = 1 : size_map(3)
        calib = squeeze(data_calib_fro(:,:,isli,:));
        [k, S] = dat2Kernel(calib, kSize);
        idx = max(find(S >= S(1)*eigThresh_1));
        [M,W] = kernelEig(k(:,:,:,1:idx), [size_map(1),size_map(2)]);
        maps = M(:,:,:,end-nset+1:end);
        weights = W(:,:,end-nset+1:end);
        weights = (weights - eigThresh_2)./(1-eigThresh_2).* (weights > eigThresh_2);
        weights = -cos(pi*weights)/2 + 1/2;
        [a, b, c] = size(weights);
        weights = repmat(reshape(weights, [a b 1 c]),[1,1,size_map(4),1]);
        smap(:,:,isli,:,:) = maps.*weights;
    end
else
    kSize3D = [7, 5, 5];
    % Compute 3D kernels
    [k,S] = dat2Kernel3D(data_calib, kSize3D);
    idx = max(find(S >= S(1)*eigThresh_1));
    [M,W] = kernelEig3D(k(:,:,:,:,1:idx), [size_map(1), size_map(2),size_map(3)]);
    
    % Compute 3D Soft-SENSE ESPIRiT Maps
    maps = M(:,:,:,:,end-nset+1:end);
    weights = W(:,:,:,end-nset+1:end);
    weights = (weights - eigThresh_2)./(1-eigThresh_2).* (weights > eigThresh_2);
    weights = -cos(pi*weights)/2 + 1/2;
    [a, b, c, d] = size(weights);
    weights = repmat(reshape(weights, [a b c 1 d]),[1,1,1,size_map(4),1]);
    smap = maps.*weights;
    %smap = M(:,:,:,:,end-nset+1:end) .* repmat(W(:,:,:,end-nset+1:end) > eigThresh_2,[1,1,1,size_map(4)]);
end

toc

end