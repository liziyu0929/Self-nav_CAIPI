function mask = samp_shot(nx,ny,nz,kysft,kzsft,p,psft,kz_cover)
mask = zeros(nx,ny,nz);
cnt = 1;
% p = [0:1:kz_cover, kz_cover-1:-1:1]; 
p = circshift(p,psft);
for ii = kysft:3:ny
    mask(:, ii, p(mod(cnt-1, 2*kz_cover)+1) + kzsft) = 1; 
    cnt = cnt + 1;
end
end