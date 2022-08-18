function X = image2patch(image, dim)
%MODULE_IM2PATCH Summary of this function goes here
% Goal : decompose the 3D tensor (e.g., video) into 2D patches
[m, n, k, nframe] = size(image);
XX = cell(dim*dim,k,nframe);
N   =   m-dim+1;
M   =   n-dim+1;
L     =   N*M;   
for nn = 1:nframe
    kk=0;
    for i  = 1:dim
        for j  = 1:dim
            kk=kk+1;
            for p = 1:k      
                blk  =  image(i : m-dim+i, j : n-dim+j, p, nn);
                XX{kk,p,nn}=blk(:);
            end
        end
    end
end

XX=reshape(XX,[],nframe);
XX=cell2mat(XX);
XX=reshape(XX,L,[],nframe);
X=permute(XX,[2 1 3]);
end


