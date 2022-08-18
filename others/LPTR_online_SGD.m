function Xhat=LPTR_online_SGD(MissM,Mask,I,option)

%% parameter setting
border_ext=option.border_ext;
border_rm=option.border_rm;
interval=option.interval;
dim=option.dim;
dims = dim/interval; 
if mod(dim,interval)~=0
    err('interval must be divided by patch size');
end
stride=option.stride;
se=option.se;
nChannel=option.nChannel;
buffersize=option.buffersize;
p_frame=option.p_frame;
batch_start_size=option.batch_start_size;

param.search_window_size=ceil(option.search_window_size/interval);
param.K=p_frame;

para_TR.r = option.TRrank*ones(ndims(MissM),1);
para_TR.max_iter = option.TRmaxiter;
para_TR.max_tot  = option.TRtot;
para_TR.disp = option.TRdisp;

[m, n, ~]=size(MissM);
border_ext_row_left=border_ext;
border_ext_row_right=border_ext-mod(n+border_ext*2,interval);
border_ext_col_up=border_ext;
border_ext_col_down=border_ext-mod(m+border_ext*2,interval);

if nChannel==1
    MissM=reshape(MissM,m,n,nChannel,[]);
    Mask=reshape(Mask,m,n,nChannel,[]);
end

%% initial
total_frame_num=size(MissM,ndims(MissM));
m=m+border_ext_col_up+border_ext_col_down;
n=n+border_ext_row_left+border_ext_row_right;
sm=kron(1:interval,ones(1,interval));
sn=kron(ones(1,interval),1:interval);
patchnew=cell(1,4);
patchmasknew=cell(1,4);

%% video patch index
total_patch_size = [m-dim+interval, n-dim+interval];
%total_buffer_size = [m-dim+interval, n-dim+interval, buffersize];
patch_ind_adjust = (m-dim+interval) * (n-dim+interval);
total_patch_size_s = [(m-dim+interval)/interval, (n-dim+interval)/interval];

%% initial patch index
m2=m-dim+interval;
n2=n-dim+interval;
border_ext2=floor(border_ext/2);

ind_ref_m =  border_ext2:stride:m2-border_ext2; 
if ind_ref_m(end)<m2-border_ext2
    ind_ref_m=[ind_ref_m m2-border_ext2];
end
ind_ref_n  =  border_ext2:stride:n2-border_ext2;
if ind_ref_n(end)<n2-border_ext2
    ind_ref_n=[ind_ref_n n2-border_ext2];
end

ind_initial=[];
for  i  =  1 : length(ind_ref_m)
    for  j  =  1 : length(ind_ref_n)
        row    =   ind_ref_m(i);
        col     =   ind_ref_n(j);
        ind_initial = [ind_initial sub2ind(total_patch_size, row, col)];
    end
end

num_patch  =  (m-dim+interval) * (n-dim+interval); 
idx_patch  =  1 : num_patch;    
idx_patch   =  reshape(idx_patch, m-dim+interval, n-dim+interval);
for i=1:interval^2   
    si=sm(i);sj=sn(i);     
    cc=mod(si-1,interval)+mod(sj-1,interval)*interval+1;
    blk_arr_temp=idx_patch(si:interval:end,sj:interval:end,:);
    blk_arr_part{cc}=blk_arr_temp(:);
end

%% initial variables
C=zeros(m,n);
blk_arr=zeros(length(ind_initial),p_frame);
Utr_TR=cell(1,length(ind_initial));

for t=1:1:total_frame_num
    
    tic
    
    if nChannel==1
        MissM_t=MissM(:,:,t);
        Mask_t=Mask(:,:,t);
    elseif nChannel==3
        MissM_t=MissM(:,:,:,t);
        Mask_t=Mask(:,:,:,t);
    end
    MissM_t = padarray(MissM_t,[border_ext_col_up border_ext_row_left],'symmetric','pre');
    MissM_t = padarray(MissM_t,[border_ext_col_down border_ext_row_right],'symmetric','post');
    Mask_t = padarray(Mask_t,[border_ext_col_up border_ext_row_left],'symmetric','pre');
    Mask_t = padarray(Mask_t,[border_ext_col_down border_ext_row_right],'symmetric','post');
    MissMse=imdilate(MissM_t,se);
    [m, n, ~]=size(MissM_t); 
    MissMdouble=double(MissM_t)/255;
    MissMdouble(Mask_t==0)=-1;
    Maskse=imdilate(Mask_t,se);
    Maskse=Maskse==1;
    
    blk_arr_pre=blk_arr;
        
    if t==1
        patch_ind = ind_initial;
        for i=1:interval^2   
            si=sm(i);sj=sn(i);
            cc=mod(si-1,interval)+mod(sj-1,interval)*interval+1;
            patchnew{cc} = image2patch(MissMse(si:interval:end,sj:interval:end,:), dims);   
            patchmasknew{cc} = image2patch(Maskse(si:interval:end,sj:interval:end,:), dims);   
            %I_s{i}=I_t(si:interval:end,sj:interval:end,:);
        end
    else
        % remove zero columns of blk_arr
        ind_zero=sum(squeeze(blk_arr(1:length(Utr_TR),1)))~=0;
        Utr_TR=Utr_TR(ind_zero);
        ind_zero=sum(blk_arr(:,1))~=0;
        ind_new=ind_new(ind_zero);
        blk_arr_pre=blk_arr_pre(ind_zero,:);
        blk_arr=blk_arr(ind_zero,:); 
        % adjust patch index and 
        blk_arr=(blk_arr-patch_ind_adjust).*(blk_arr~=0);
        patch_ind = ind_new;
    end
    
    C_pre=C;
    patchold = patchnew;
    patchmaskold = patchmasknew;
    ind_new=zeros(1,length(patch_ind));
    patch_ind_s=zeros(1,length(patch_ind));
    patch_ind_c=zeros(1,length(patch_ind));
       
    %% compute new patch locations
    if t>1
        for i=1:interval^2   
            si=sm(i);sj=sn(i);
            cc=mod(si-1,interval)+mod(sj-1,interval)*interval+1;
            patchnew{cc} =  image2patch(MissMse(si:interval:end,sj:interval:end,:), dims);
            patchmasknew{cc} = image2patch(Maskse(si:interval:end,sj:interval:end,:), dims);
        end
    end
    
    for i =1:length(patch_ind)        
        if patch_ind(i)~=0
            [xx,yy]=ind2sub(total_patch_size,patch_ind(i));
            xxn=floor((xx-1)/interval)+1;
            yyn=floor((yy-1)/interval)+1;
            patch_ind_c(i)=mod(xx-1,interval)+mod(yy-1,interval)*interval+1;
            patch_ind_s(i)=sub2ind(total_patch_size_s, xxn, yyn);           
        end    
    end
    
    blk_patch=zeros(interval^2*p_frame,length(patch_ind_s));
    err_patch=blk_patch;
    
    for i =1:length(patch_ind_s)        
       if patch_ind_s(i)~=0
           refpatch = patchold{patch_ind_c(i)}(:,patch_ind_s(i));
           refmask = patchmaskold{patch_ind_c(i)}(:,patch_ind_s(i));
           for s=1:1:interval^2                 
                [blk_patch_s,error_arr_s] = patch_BM_adjacent_frame(patchnew{s}, patchmasknew{s}, refpatch, refmask, total_patch_size_s, param, patch_ind_s(i)); 
                blk_patch((s-1)*p_frame+1:s*p_frame,i)=blk_arr_part{s}(blk_patch_s);
                err_patch((s-1)*p_frame+1:s*p_frame,i)=error_arr_s;
           end
       end
    end
    
    for i=1:length(patch_ind)        
        if patch_ind(i)~=0
            [xx,yy]=ind2sub(total_patch_size,patch_ind(i));
            xxn=xx:min(xx+dim-1,m);
            yyn=yy:min(yy+dim-1,n);
            if t>5&&min(min(C_pre(xxn,yyn)))>6&&sum(squeeze(blk_arr_pre(:,i,1))~=0)==buffersize
               continue;
            end
            [~, ind] = sort(err_patch(:,i));
            if t==1&&blk_patch(ind(1),i)~=patch_ind(i)
                blk_arr(i,1)=patch_ind(i);
                blk_arr(i,2:p_frame)=blk_patch(ind(1:p_frame-1),i);
                ind_new(i)=patch_ind(i);
            else
                blk_arr(i,:)=blk_patch(ind(1:p_frame),i);
                ind_new(i)=blk_patch(ind(1),i);
            end           
       end
    end
    
    %% compute overlaps
    C=zeros(m,n);
    for i =1:length(ind_new)      
        if ind_new(i)~=0      
            [xx,yy]=ind2sub(total_patch_size,ind_new(i));
            xxn=xx:min(xx+dim-1,m);
            yyn=yy:min(yy+dim-1,n);        
            C(xxn,yyn)=C(xxn,yyn)+ones(length(xxn),length(yyn));              
        end
    end
          
    %% detect holes and add patches
    BW=imbinarize(C);
    [L,num]= bwlabel(~BW,8);   
    stats = regionprops(L,'BoundingBox');
    
    ind_add=[];
    
    for j=1:1:num
        bbx=struct2array(stats(j));
        ind_ref_m  =  max(floor(bbx(2)-border_ext2),border_ext2):stride:min(floor(bbx(2)+bbx(4)+border_ext2),m2-border_ext2); 
        ind_ref_n   =  max(floor(bbx(1)-border_ext2),border_ext2):stride:min(floor(bbx(1)+bbx(3)+border_ext2),n2-border_ext2);
        if ~isempty(ind_ref_m)&&~isempty(ind_ref_n)
            if ind_ref_m(end)<m2-border_ext2
                ind_ref_m=[ind_ref_m m2-border_ext2];
            end
            if ind_ref_n(end)<n2-border_ext2
                ind_ref_n=[ind_ref_n n2-border_ext2];
            end
        end
        for  ii  =  1 : length(ind_ref_m)
            for  jj  =  1 : length(ind_ref_n)
            row    =   ind_ref_m(ii);
            col     =   ind_ref_n(jj);
            patchbw=BW(row:min(row+dim-1,m),col:min(col+dim-1,n));
            if sum(patchbw(:))<dim^2
                ind_add = [ind_add sub2ind(total_patch_size, row, col)];
            end
            end
        end
    end

    ind_add_s=zeros(1,length(ind_add));
    ind_add_c=zeros(1,length(ind_add));
    
    for i =1:length(ind_add)        
        if ind_add(i)~=0
            [xx,yy]=ind2sub(total_patch_size,ind_add(i));
            xxn=floor((xx-1)/interval)+1;
            yyn=floor((yy-1)/interval)+1;
            ind_add_c(i)=mod(xx-1,interval)+mod(yy-1,interval)*interval+1;
            ind_add_s(i)=sub2ind(total_patch_size_s, xxn, yyn); 
        end
    end
    
    blk_patch=zeros(interval^2*p_frame,length(ind_add_s));
    err_patch=blk_patch;
    
    for i =1:length(ind_add_s)        
       refpatch = patchnew{ind_add_c(i)}(:,ind_add_s(i));
       refmask = patchmasknew{ind_add_c(i)}(:,ind_add_s(i));
       for s=1:1:interval^2                 
            [blk_patch_s,error_arr_s] = patch_BM_adjacent_frame(patchnew{s}, patchmasknew{s}, refpatch, refmask, total_patch_size_s, param, ind_add_s(i)); 
            blk_patch((s-1)*p_frame+1:s*p_frame,i)=blk_arr_part{s}(blk_patch_s);
            err_patch((s-1)*p_frame+1:s*p_frame,i)=error_arr_s;
       end
    end
    
    for i=1:length(ind_add)        
        [xx,yy]=ind2sub(total_patch_size,ind_add(i));
        xxn=xx:min(xx+dim-1,m);
        yyn=yy:min(yy+dim-1,n);
        C(xxn,yyn)=C(xxn,yyn)+ones(length(xxn),length(yyn)); 
        [~, ind] = sort(err_patch(:,i));
        if blk_patch(ind(1),i)~=ind_add(i)
            blk_arr(i+length(patch_ind),1)=ind_add(i);
            blk_arr(i+length(patch_ind),2:p_frame)=blk_patch(ind(1:p_frame-1),i);
        else
            blk_arr(i+length(patch_ind),:)=blk_patch(ind(1:p_frame),i);  
        end    
        
    end
          
    ind_new=[ind_new ind_add];
     
    Utr_TR_pre=Utr_TR;
    Utr_TR=cell(1,length(ind_new));
    Utr_TR(1:length(Utr_TR_pre))=Utr_TR_pre;

    %% update all subtensors  
    for kk=1:length(ind_new)
        if ind_new(kk)~=0
            arr=squeeze(blk_arr(kk,:))';
            arr=arr(:);               
            n_p=min(p_frame,sum(arr~=0));
            if nChannel==1
                patch_miss=zeros(dim,dim,n_p);   
            else
                patch_miss=zeros(dim,dim,nChannel,n_p);  
            end
            p_count=1;
            for i=1:1:length(arr)
                if arr(i)==0
                    continue;
                end
                [xx,yy]=ind2sub(total_patch_size,arr(i));
                xxn=xx:min(xx+dim-1,m);
                yyn=yy:min(yy+dim-1,n);
                if nChannel==1
                    patch_miss(1:length(xxn),1:length(yyn),p_count)=MissMdouble(xxn,yyn);
                elseif nChannel==3
                    patch_miss(1:length(xxn),1:length(yyn),:,p_count)=MissMdouble(xxn,yyn,:);
                end
                p_count=p_count+1;
            end
            if ~isempty(Utr_TR{kk})
                Utr_TR{kk} = TRSGD_scaled_add_acc(patch_miss, double(patch_miss~=-1), para_TR, Utr_TR{kk}, p_frame); 
            else 
                Utr_TR{kk} = TRSGD_scaled_initial(patch_miss, double(patch_miss~=-1), para_TR);        
            end
        end
    end

    
    %% output the first image in the buffer         
    I_hat=zeros(m,n,nChannel);
    C_hat=zeros(m,n);

    blk_arr_out=blk_arr(:,1);

    blk_arr_out_ind=find(blk_arr_out~=0);

    for kk=1:1:length(blk_arr_out_ind)

        Utr_n=Utr_TR{blk_arr_out_ind(kk)};         
        patch_hat=reshape(fullTR(Utr_n),dim,dim,nChannel,p_frame);  
        for ii=1:1:p_frame
            [xx,yy]=ind2sub(total_patch_size,blk_arr(blk_arr_out_ind(kk),ii));
            xxn=xx:min(xx+dim-1,m);
            yyn=yy:min(yy+dim-1,n);
            I_hat(xxn,yyn,:)=I_hat(xxn,yyn,:)+patch_hat(1:length(xxn),1:length(yyn),:,ii);   
            C_hat(xxn,yyn)=C_hat(xxn,yyn)+ones(length(xxn),length(yyn));  
        end
     end
%         
%         
%         if t-buffersize+1==1
%             fprintf('Output frame, time cost, number of patches: ');
%         else
%             fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
%         end
%         %fprintf(['Output frame: ', num2str(t-buffersize+1), ', time cost: ' num2str(toc) 's', ', number of patches: ' num2str(length(blk_arr_out_ind))]);    
%         fprintf('%3i, %.4fs, %5i', t-buffersize+1, toc, length(blk_arr_out_ind));     
%         %imshow(I_hat(border_ext+1:end-border_ext,border_ext+1:end-border_ext,:)./C_hat(border_ext+1:end-border_ext,border_ext+1:end-border_ext))
%         %Xhat(:,:,:,t-buffersize+1)=I_hat(border_ext+1:end-border_ext,border_ext+1:end-border_ext,:)./C_hat(border_ext+1:end-border_ext,border_ext+1:end-border_ext);
    Xhat(:,:,:,t)=I_hat(border_ext_col_up+1:end-border_ext_col_down,border_ext_row_left+1:end-border_ext_row_right,:)./C_hat(border_ext_col_up+1:end-border_ext_col_down,border_ext_row_left+1:end-border_ext_row_right);
    disp(psnr(I(:,:,:,t),Xhat(:,:,:,t)));
        

                
end

fprintf('\n');

Xhat(isnan(Xhat))=0;

Xhat=squeeze(Xhat);

end

