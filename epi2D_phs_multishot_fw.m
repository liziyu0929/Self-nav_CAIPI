% forward model for multi-shot EPI reconstruction
%
% Ziyu Li, Wenchuan Wu, University of Oxford, 2024

function dataout = epi2D_phs_multishot_fw(imgin,phs,sampling_order,z_start)

% assumption: phase is constant within slab for each shot

[Ny,Nz,Nch] = size(imgin);

n_shot = size(sampling_order,3);
dataout = zeros(Ny,Nz,Nch,n_shot);


for i_shot = 1 : n_shot

    phs_shot = squeeze(phs(:,i_shot,:));
    phs_shot = repmat(reshape(phs_shot,[Ny, 1, Nch]),[1,Nz,1]);
    
    datatmp = fft2c(bsxfun(@times,imgin,phs_shot));
    ky_curr = sampling_order(:,1,i_shot);
    kz_curr = sampling_order(:,2,i_shot);
    
    indy = ky_curr(find(ky_curr(:)));
    indz = kz_curr(find(kz_curr(:)))-z_start+1;
    
    for iii = 1 : length(indy)
        dataout(indy(iii),indz(iii),:,i_shot) =  datatmp(indy(iii),indz(iii),:);
    end

end
    
end

