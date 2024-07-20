function mask = samp2mask(ny, nz, sampling_order, disp_flag)
mask2show = zeros(ny,nz);
if disp_flag
    figure;
end
mask = zeros(ny,nz);
for i_ex = 1 : size(sampling_order, 3)
    for iy = 1 : size(sampling_order, 1)
        if sampling_order(iy,1,i_ex) * sampling_order(iy,2,i_ex) > 0
            mask2show(sampling_order(iy,1,i_ex),sampling_order(iy,2,i_ex)) = 1;  
            mask(sampling_order(iy,1,i_ex),sampling_order(iy,2,i_ex)) = 1;  
        end
        if disp_flag
            imshow(mask2show)
            title(['shot number = ',num2str(i_ex)])
            pause(0.05);
        end
    end
end
