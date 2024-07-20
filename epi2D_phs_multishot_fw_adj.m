% adjoint of forward model for multi-shot EPI reconstruction
%
% Ziyu Li, Wenchuan Wu, University of Oxford, 2024

function dataout = epi2D_phs_multishot_fw_adj(kin,phs,sampling_order,z_start)

% assumption: phase is constant within slab for each shot

[Ny,Nz,Nch,n_shot] = size(kin);
dataout = zeros(Ny,Nz,Nch);


for i_shot = 1 : n_shot
    
    phs_shot = squeeze(phs(:,i_shot,:));
    phs_shot = repmat(reshape(phs_shot,[Ny, 1, Nch]),[1,Nz,1]);
    
    data_tmp = zeros(Ny,Nz,Nch);
    
    ky_curr = sampling_order(:,1,i_shot);
    kz_curr = sampling_order(:,2,i_shot);
    indy = ky_curr(find(ky_curr(:)));
    indz = kz_curr(find(kz_curr(:)))-z_start+1;
    
    for iii = 1 : length(indy)
   
        data_tmp(indy(iii),indz(iii),:) = kin(indy(iii),indz(iii),:,i_shot);
        % kin(indy(iii),indz(iii),:) = 0;
    end
    % data_tmp(indy,indz,:)=kin(indy,indz,:);
    
    im_tmp=bsxfun(@times,ifft2c(data_tmp),conj(phs_shot));

    dataout = dataout + im_tmp;

end

    
end

