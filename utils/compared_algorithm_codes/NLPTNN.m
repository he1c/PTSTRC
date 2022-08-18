function Xhat=NLPTNN(MissM,Mask,I,option)

radius = option.radius;
search_rad = option.search_rad;
search_gap = option.search_gap;
THRESHOLD = option.THRESHOLD;
nChannel=option.nChannel;

[m, n, ~]=size(MissM);
if nChannel==1
    MissM=reshape(MissM,m,n,nChannel,[]);
    Mask=reshape(Mask,m,n,nChannel,[]);
end
num_frame=size(MissM,4);
Xhat = zeros(size(MissM)); 

if num_frame>=3
    MissM=cat(4,MissM(:,:,:,3),MissM,MissM(:,:,:,end-2));
    Mask=cat(4,Mask(:,:,:,3),Mask,Mask(:,:,:,end-2));
    fprintf('NLPTNN Frame, time cost:     ');
    for k=2:1:num_frame+1
        fprintf('\b\b\b%3i',k-1);
        for i = 1:nChannel
            Xhat(:,:,i,k-1) = inpaint_multi(MissM(:,:,i,k-1:k+1), ~Mask(:,:,i,k-1:k+1), radius, search_rad, search_gap, @t_SVD_121norm, THRESHOLD);
        end
        if exist(I)
            disp(psnr(I(:,:,:,k-1),(1-Mask(:,:,:,k-1)).*Xhat(:,:,:,k-1)+Mask(:,:,:,k-1).*MissM(:,:,:,k-1)));
        end
    end
    fprintf('\n');
elseif num_frame<3
    for k=1:1:num_frame
        for i = 1:nChannel
            Xhat(:,:,i,k) = inpaint_single(MissM(:,:,i,k), ~Mask(:,:,i,k), radius, search_rad, search_gap, @t_SVD_121norm, THRESHOLD);
        end
    end
end

Xhat=squeeze(Xhat);

end


