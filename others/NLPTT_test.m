function Xhat=NLPTT_test(MissM,Mask,I,option)

K=option.K;
dim=option.dim;
stride=option.stride; 
nChannel=option.nChannel;
%n_patch=option.n_patch;

para_TT.beta=option.beta;
para_TT.maxitr=option.maxitr;
para_TT.stopc=option.stopc;
para_TT.debug=0;

%% NLPTT
MissM2=MissM;
MissM2(Mask==0)=-1;
Mask=logical(Mask);
sz=size(MissM);
m=sz(1);n=sz(2);
framenum=option.framenum;
MissM_int=uint8(MissM*255);

ind_ref_m =  20:stride:m-dim+1; 
if ind_ref_m(end)<m-dim+1
    ind_ref_m=[ind_ref_m m-dim+1];
end
ind_ref_n  =  20:stride:n-dim+1;
if ind_ref_n(end)<n-dim+1
    ind_ref_n=[ind_ref_n n-dim+1];
end

zz=length(ind_ref_m)*length(ind_ref_n);

patchrefall=uint8(zeros(dim^2*prod(sz(3:end)),zz));
patchrefMaskall=logical(patchrefall);

for  i  =  1 : length(ind_ref_m)
    for  j  =  1 : length(ind_ref_n)
        row    =   ind_ref_m(i);
        col     =   ind_ref_n(j);
        patchref  =  MissM_int(row:row+dim-1, col:col+dim-1,:);
        patchrefMask  =  Mask(row:row+dim-1, col:col+dim-1,:);
        patchrefall(:,(i-1)*length(ind_ref_n)+j)=patchref(:);
        patchrefMaskall(:,(i-1)*length(ind_ref_n)+j)=patchrefMask(:);
    end
end


metric=zeros((m-dim+1)*(n-dim+1),zz);

pt_ind=[];

for i=1:m-dim+1
    for j=1:n-dim+1
        pt_ind=[pt_ind sub2ind([m-dim+1,n-dim+1],i,j)];
    end
end

metric2=zeros(length(pt_ind),zz);
metric_ind=zeros(length(pt_ind),1);

tic

parfor mm=1:length(pt_ind)
    [row,col]=ind2sub([m-dim+1,n-dim+1],pt_ind(mm));
    patchcand =  MissM_int(row:row+dim-1, col:col+dim-1,:);
    patchcandMask  =  Mask(row:row+dim-1, col:col+dim-1,:);
    patchcandMask = patchcandMask(:);
    obind=patchcandMask==1;
    patchcand=patchcand(:);
    patchcand=patchcand(obind);
    patchcandMask=patchcandMask(obind);
    patchrefMask=patchrefMaskall(obind,:);
    patchref=patchrefall(obind,:);
    for kk=1:1:zz
        dis2=single(patchcand)-single(patchref(:,kk));
        coMask=single(patchcandMask.*patchrefMask(:,kk));
        metric2(mm,kk)=(sum((coMask.*dis2).^2)+1e-5)/(sum(coMask)+1e-10);
        metric_ind(mm)=pt_ind(mm);
    end
end

toc

disp('PM end');

metric(metric_ind,:)=metric2;

[error, ind]=sort(metric);

% [x,y]=meshgrid(ind_ref_m,ind_ref_n);
% indall=(sub2ind([m-dim+1,n-dim+1],x,y));
% indall=indall(:);
% sum(abs(ind(1,:)-indall))

blk_arr_n=ind(1:K,:);

%% aggreagate
I_hat=zeros(m,n,nChannel,framenum);
C_hat=zeros(m,n,1);
patch_hat=cell(1,size(blk_arr_n,2));

parfor kk=1:1:size(blk_arr_n,2)
    disp(kk)
    arr=blk_arr_n(:,kk);
    patch_miss=zeros(dim,dim,nChannel,framenum,K);
    patch_full=zeros(dim,dim,nChannel,framenum,K);
    for i=1:1:K
        [xx,yy]=ind2sub([m-dim+1,n-dim+1],arr(i));
        xxn=xx:min(xx+dim-1,m);
        yyn=yy:min(yy+dim-1,n);
        patch_miss(1:length(xxn),1:length(yyn),:,:,i)=MissM2(xxn,yyn,:,:); 
        patch_full(1:length(xxn),1:length(yyn),:,:,i)=I(xxn,yyn,:,:); 
    end
    %patch_miss_t=reshape(patch_miss,n_patch);
    patch_miss_t=squeeze(patch_miss);

    tic
    patch_rec=TTNN(patch_miss_t, double(patch_miss_t~=-1), [], para_TT);
    patch_hat{kk}=reshape(patch_rec,dim,dim,nChannel,framenum,K);
    toc
    
    tic
    patch_rec=TTNN(patch_miss_t(:,:,:,:,1:10), double(patch_miss_t(:,:,:,:,1:10)~=-1), [], para_TT);
    patch_hat2=reshape(patch_rec,dim,dim,nChannel,framenum,10);
    toc
    
    disp([psnr(patch_hat{kk},patch_full) psnr(patch_hat2,patch_full(:,:,:,:,1:10))])
    

end

% for kk=1:1:size(blk_arr_n,2)
% 
%       for i=1:1:K
%           [xx,yy]=ind2sub([m-dim+1,n-dim+1],blk_arr_n(i,kk));
%           xxn=xx:min(xx+dim-1,m);
%           yyn=yy:min(yy+dim-1,n);
%           I_hat(xxn,yyn,:,:)=I_hat(xxn,yyn,:,:)+patch_hat{kk}(1:length(xxn),1:length(yyn),:,:,i);   
%           C_hat(xxn,yyn)=C_hat(xxn,yyn)+ones(length(xxn),length(yyn));  
%       end
% 
% end
% 
Xhat=0;%squeeze(I_hat./C_hat);

end



        