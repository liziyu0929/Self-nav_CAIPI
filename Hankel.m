function [H] =Hankel(k,r)
%   Hankel matrix construction

 % r: filter size

if ndims(k)==3
    k = permute(k,[1,2,4,3]);% m*n*coil*shot
end
[m,n,coil,shot]=size(k);
H=zeros((m-r+1)*(n-r+1), r*r*coil, shot);
for s=1:shot
    for c=1:coil
        for j=1:n-r+1
            for i=1:m-r+1
                  H(i+(j-1)*(m-r+1),r*r*(c-1)+1:r*r*c,s)=...
                      reshape(k(i:i+r-1,j:j+r-1,c,s),1,r*r);
            end
        end
    end
end
 
H = reshape(H, [(m-r+1)*(n-r+1), r*r*coil*shot]);
 
end