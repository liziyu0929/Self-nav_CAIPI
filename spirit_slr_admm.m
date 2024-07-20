% SPIRiT structured low-rank reconstruction of 2D phase map using ADMM
%
% Ziyu Li, Wenchuan Wu, University of Oxford, 2024

function res = spirit_slr_admm(y, niter, samp, GOP, x0, x_mag)

% set parameters
[sx,sy,nc,nshot] = size(y);
r = 10; % Hankel matrix kernel size
lambda1 = 1; lambda2 = nshot*nc; alpha = 1e-4;
rho = 1e-4;

x = x0;
x_mag = repmat(x_mag,[1,1,1,nshot]);
x1 = fft2c(x_mag .* exp(1i.*angle(ifft2c(x0))));
z = Hankel(x, r);
Hx = z;
u = z*0;

[~, N]  =   pinv_hankel(Hankel(ones(sx,sy,nc,nshot), r), r, sx, sy, nc, nshot);
N = squeeze(N);

H_fwd   =   @(x)Hankel(x, r);
H_adj   =   @(x)squeeze(pinv_hankel(x, r, sx, sy, nc, nshot)).*N;

% argmin_x {||Ax-y||_2^2 + \lambda_1||(G-I)x||_2^2 + \lambda_2||z||_1 +
% alpha||x - x1||_2^2}
% s.t. z = Hx

for i = 1:niter
    disp(['iter:' num2str(i)])
    % z subproblem
    % argmin_z { lambda||z||_* + (rho/2)||z - Hx + u||_2^2}
    % singular value thresholding
    tic
    [U,S,V] = svd(single(Hx - u), 'econ');
    % [U,S,V] = svd(Hx - u, 'econ');
    keep = 1 : size(S, 1)/lambda2;
    z = U(:,keep) * S(keep,keep) * V(:,keep)';
    clear U S V
    toc
    
    % x subproblem
    % argmin_x { (1/2)||Ax-y||_2^2 + (lambda_1/2)||(G-I)x||_2^2 
    % + (rho/2)||z - Hx + u||_2^2 + (alpha/2)||x - x1||_2^2}
    % (A'A + (lambda_1)(G-I)'(G-I) + (rho)H'H + (alpha)I)x = A'y + (rho)H'(z+u) + alpha*x1
    tic
    y_tmp = reshape((1/2)*A(y, samp)+ (rho/2)*H_adj(z+u) + alpha*x1 ,[],1);
    [tmp,~]   =   pcg(@(x)cgfun(x, sx, sy, nc, samp, GOP, lambda1, rho, alpha, N),...
                    y_tmp, 1E-5, 30, [], [], x0(:));
    x       =   reshape(tmp, sx, sy, nc, nshot);
    clear tmp x0

    Hx = H_fwd(x);
    toc
    
    % dual update
    u       =   u + z - Hx;
    x0 = x;
    x1 = fft2c(x_mag .* exp(1i.*angle(ifft2c(x0))));
    
end

res  = ifft2c(reshape(x,[sx,sy,nc,nshot]));

end

% A_fwd = A_adj
function res = A(x, samp)
res  = x.*samp;
end

function q = cgfun(x, sx, sy, nc, samp, GOP, lambda, rho, alpha, N)

x = reshape(x, sx, sy, nc, []); nshot = size(x, 4);
spirit_reg = zeros(sx, sy, nc, nshot);
for ishot = 1:nshot
    x_shot = squeeze(x(:,:,:,ishot));
    samp_shot = squeeze(samp(:,:,:,ishot));
    spirit_reg(:,:,:,ishot) = (1/2)*A(A(x_shot, samp_shot), samp_shot) + ...
        (lambda/2) .* (GOP'*(GOP*x_shot));
end
q = spirit_reg + (rho/2)*(N.*x) + alpha*x; 
q = reshape(q,[],1);
end