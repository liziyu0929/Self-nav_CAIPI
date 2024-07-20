function [res, RESVEC] = cgL1SPIRiT(y,x0, GOP, nIterCG, XOP, lambda, alpha,nIterSplit)
%
%
%  res = cgSPIRiT(y,GOP, nIter, lambda,x0)
%  
%  Implementation of the Cartesian, conjugate gradiend SPIRiT reconstruction
%
%  Input:
%		y -	Undersampled k-space data. Make sure that empty entries are zero
%			or else they will not be filled.
%		GOP -	the SPIRiT kernel operator obtained by calibration
%		nIter -	Maximum number of iterations
%		lambda-	Tykhonov regularization parameter
%       x0    - Initil value
%
% Outputs:
%		res - Full k-space matrix
%
%
%

N = numel(x0);
M = numel(y);
imSize = size(x0);

% make dyadic size if Wavelets are used. 
if strcmp(class(XOP),'Wavelet') == 1

    if length(imSize)>2
        imSize_dyd = [max(2.^ceil(log2(imSize(1:2)))), max(2.^ceil(log2(imSize(1:2)))),imSize(3)];
    else
        imSize_dyd = [max(2.^ceil(log2(imSize(1:2)))), max(2.^ceil(log2(imSize(1:2))))];
    end
else
    imSize_dyd = imSize;
end

kernel = getKernel(GOP);
kSize = [size(kernel,1),size(kernel,2)];

[sx,sy,nCoils] = size(y);

idx_acq = find(abs(y)>0);
idx_nacq = find(abs(y)==0);
N = length(idx_nacq(:));


yy = GOP*y; 
dataSize = [size(y)];

res = x0(:);

for n=1:nIterSplit;
    b = [-yy(:); sqrt(alpha)*res(idx_nacq)];

    [tmpres,FLAG,RELRES,ITER,RESVEC] = lsqr(@aprod,b,1e-6,nIterCG, speye(N,N),speye(N,N),x0(idx_nacq),GOP,sx,sy,nCoils,idx_nacq, alpha);
    res = y;
    res(idx_nacq) = tmpres;
    
    if lambda > 0
        tmp = ifft2c(res);
        tmp = zpad(tmp,imSize_dyd);
        tmp = XOP*tmp;
        tmp = SoftThresh(tmp,lambda/sqrt(alpha));
        res = XOP'*tmp;
        res = reshape(res,imSize_dyd);
        res = crop(res,imSize);
        res = fft2c(res);
    end
                                                                
end




function [res,tflag] = aprod(x,GOP,sx,sy,nCoils,idx_nacq, alpha,tflag)
	
	kernel = getKernel(GOP);
	kSize = [size(kernel,1),size(kernel,2)];

	if strcmp(tflag,'transp');
		tmpy = reshape(x(1:sx*sy*nCoils),sx,sy,nCoils);
        res = GOP'*tmpy;
        res = res(idx_nacq)+ x(sx*sy*nCoils+1:end)*sqrt(alpha);
    else
        tmpx = zeros(sx,sy,nCoils);
		tmpx(idx_nacq) = x;
		res = GOP*tmpx;
		res = [res(:) ; sqrt(alpha)*x(:)];
	end






