function Xhat=LPTR_mix_test(MissM,Mask,I,option)

%% parameter setting
border_ext=option.border_ext;
interval=option.interval;
dim=option.dim;
K=option.K;
K_TRNN=option.K_TRNN;
thre_TRNN=option.thre_TRNN;
param.search_window_size=ceil(option.search_window_size/interval);
dims  = dim/interval;
if mod(dim,interval)~=0
    err('interval must be divided by patch size');
end
param.K=K;
param.stride=floor(option.stride/interval); 
param.border_ext=floor(border_ext/interval);
se=option.se;
nChannel=option.nChannel;
adjnum =option.adjnum;
buffersize=2*adjnum+1;
refind=adjnum+1;

para_TR.r = option.TRrank*ones(4,1);
para_TR.max_iter = option.TRmaxiter;
para_TR.max_tot  = option.TRtot;
para_TR.disp = option.TRdisp;

para_TRNN.beta=option.beta;
para_TRNN.alpha=option.alpha*ones(4,1);
para_TRNN.stopc=option.stopc;
para_TRNN.maxitr=option.maxitr;
para_TRNN.debug=option.debug;

%% extend video sequence
[m, n, ~]=size(MissM);
if nChannel==1
    MissM=reshape(MissM,m,n,nChannel,[]);
    Mask=reshape(Mask,m,n,nChannel,[]);
    para_TR.r = option.TRrank*ones(3,1);
end
total_frame_num=size(MissM,4);
Xhat = double(MissM)/255;
MissM=cat(4,MissM(:,:,:,adjnum+1:-1:2),MissM,MissM(:,:,:,end-adjnum:end-1));
Mask=cat(4,Mask(:,:,:,adjnum+1:-1:2),Mask,Mask(:,:,:,end-adjnum:end-1));

%% pre-processing
border_ext_row_left=border_ext;
border_ext_row_right=border_ext-mod(n+border_ext*2,interval);
border_ext_col_up=border_ext;
border_ext_col_down=border_ext-mod(m+border_ext*2,interval);
MissM = padarray(MissM,[border_ext_col_up border_ext_row_left],'symmetric','pre');
MissM = padarray(MissM,[border_ext_col_down border_ext_row_right],'symmetric','post');
Mask = padarray(Mask,[border_ext_col_up border_ext_row_left],'symmetric','pre');
Mask = padarray(Mask,[border_ext_col_down border_ext_row_right],'symmetric','post');
I = padarray(I,[border_ext_col_up border_ext_row_left],'symmetric','pre');
I = padarray(I,[border_ext_col_down border_ext_row_right],'symmetric','post');

MissMse=imdilate(MissM,se);
MissMdouble=double(MissM)/255;
Idouble=double(I)/255;
MissMdouble(Mask==0)=-1;
Maskse=imdilate(Mask,se);
Maskse=Maskse==1;
[m, n, ~]=size(MissM); 

%% initial 
sm=kron(1:interval,ones(1,interval));
sn=kron(ones(1,interval),1:interval);
blk_arr_part=cell(1,interval^2);
patchMiss=cell(1,interval^2);
patchMask=cell(1,interval^2);
    
numTotalPatch_all  =  (m-dim+interval) * (n-dim+interval) * buffersize; 
idxTotalPatch  =  1 : numTotalPatch_all;    
idxTotalPatch   =  reshape(idxTotalPatch, m-dim+interval, n-dim+interval, buffersize);

idxTotalPatch_s  =  1 : numTotalPatch_all/interval^2;    
idxTotalPatch_s   =  reshape(idxTotalPatch_s, (m-dim+interval)/interval, (n-dim+interval)/interval, buffersize);

for i=1:interval^2   
    
    si=sm(i);sj=sn(i);
    blk_arr_temp=idxTotalPatch(si:interval:end,sj:interval:end,:);
    blk_arr_part{i}=blk_arr_temp(:);
    videoMiss=MissMse(si:interval:end,sj:interval:end,:,1:buffersize);    
    videoMask=Maskse(si:interval:end,sj:interval:end,:,1:buffersize);
    patchMiss{i} = video2patch(videoMiss, dims);  
    patchMask{i} = video2patch(videoMask, dims);  
    
end

for framenum=1+adjnum:total_frame_num+adjnum
    
    tic
    
    adjind=framenum-adjnum:framenum+adjnum;
    startind=adjind(1);
    endind=adjind(end);

    error_all=[];
    blk_all=[];
    
    tic
    
    if framenum>1+adjnum
        for i=1:interval^2   
            si=sm(i);sj=sn(i);
            patchMiss_add = image2patch(MissMse(si:interval:end,sj:interval:end,:,endind), dims);     
            patchMask_add = image2patch(Maskse(si:interval:end,sj:interval:end,:,endind), dims); 
            patchMiss{i} = cat(3, patchMiss{i}, patchMiss_add);  
            patchMask{i} = cat(3, patchMask{i}, patchMask_add);  
            patchMiss{i}(:,:,1)=[];
            patchMask{i}(:,:,1)=[];
        end
    end

    patchref=patchMiss{1};
    patchMaskref=patchMask{1};

    for i=1:interval^2     
       [blk_arr,error_arr] = patch_BM_multi_frame(patchMiss{i}, patchMask{i}, patchref, patchMaskref, idxTotalPatch_s, param, refind);
       error_all=[error_all;error_arr];
       blk_all=[blk_all;reshape(blk_arr_part{i}(blk_arr(:)),size(error_arr))];
    end

    blk_arr_n=[];
    
    patch_arr_s = get_patch_ind(idxTotalPatch_s, param, refind);
    patch_arr = blk_arr_part{1}(patch_arr_s);
    
    flag_TRNN=zeros(size(blk_all,2),1);

    for i=1:1:size(blk_all,2)
        [error, ind] = sort(error_all(:,i));
        if blk_all(ind(1),i)~=patch_arr(i)
            blk_arr_n(1,i)=patch_arr(i);
            blk_arr_n(2:K,i)=blk_all(ind(1:K-1),i);
            error_sum=sum(error(1:K-1));
        else
            blk_arr_n(:,i)=blk_all(ind(1:K),i);
            error_sum=sum(error(2:K));
        end
        if error_sum/(K-1)<thre_TRNN
            flag_TRNN(i)=1;
        end
    end

    %% aggreagate
    I_hat=zeros(m,n,nChannel);
    C_hat=zeros(m,n,1);
    %Utr_TR=cell(1,size(blk_arr_n,2));
    K_n=K*ones(size(blk_arr_n,2),1);

    for kk=1:size(blk_arr_n,2)

        arr=blk_arr_n(:,kk);
        patch_miss=zeros(dim,dim,nChannel,K);  
        for i=1:1:K
            [xx,yy,zz]=ind2sub(size(idxTotalPatch),arr(i));
            xxn=xx:min(xx+dim-1,m);
            yyn=yy:min(yy+dim-1,n);
            patch_miss(1:length(xxn),1:length(yyn),:,i)=MissMdouble(xxn,yyn,:,zz+startind-1); 
            patch_full(1:length(xxn),1:length(yyn),:,i)=Idouble(xxn,yyn,:,zz+startind-1); 
        end
        patch_miss=squeeze(patch_miss);
        if flag_TRNN(kk)
            patch_miss=patch_miss(:,:,:,1:K_TRNN);
            patch_full=patch_full(:,:,:,1:K_TRNN);
            patch_hat{kk}=TRNN(patch_miss, double(patch_miss~=-1), [], para_TRNN);
            K_n(kk)=K_TRNN;
        else            
            Utr_TR = TRSGD_scaled(patch_miss, double(patch_miss~=-1), para_TR);%Completion_TR_norm(patch_miss, double(patch_miss~=-1), para_TR);
            patch_hat{kk}=reshape(fullTR(Utr_TR'),dim,dim,nChannel,K);  
        end

    end

    for kk=1:1:size(blk_arr_n,2)
          for ii=1:1:K_n(kk)
            [xx,yy,zz]=ind2sub(size(idxTotalPatch),blk_arr_n(ii,kk));
            if zz==refind
                xxn=xx:min(xx+dim-1,m);
                yyn=yy:min(yy+dim-1,n);
                I_hat(xxn,yyn,:)=I_hat(xxn,yyn,:)+patch_hat{kk}(1:length(xxn),1:length(yyn),:,ii);   
                C_hat(xxn,yyn)=C_hat(xxn,yyn)+ones(length(xxn),length(yyn));  
            end
          end

    end
    
    if framenum-adjnum==1
        fprintf('Output frame, time cost, rate of using TRNN: ');
    else
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
    end
    fprintf('%3i, %.4fs, %.4f', framenum-adjnum, toc, sum(K_n==K_TRNN)/length(K_n)); 
    %disp(['Output frame: ', num2str(framenum-adjnum), ', time cost: ' num2str(toc) 's']);
    %imshow(I_hat(border_ext+1:end-border_ext,border_ext+1:end-border_ext,:)./C_hat(border_ext+1:end-border_ext,border_ext+1:end-border_ext))
    Xhat(:,:,:,framenum-adjnum)=I_hat(border_ext_col_up+1:end-border_ext_col_down,border_ext_row_left+1:end-border_ext_row_right,:)./C_hat(border_ext_col_up+1:end-border_ext_col_down,border_ext_row_left+1:end-border_ext_row_right);
    
end

fprintf('\n');

Xhat(isnan(Xhat))=0;

Xhat=squeeze(Xhat);

end
