clear

addpath(genpath(pwd))

opt.framenum=1;

load results\image_p_new.mat

I_OPTRC_all_old=I_OPTRC_all;
time_OPTRC_old=time_OPTRC;

I_OPTRC_all=[];
time_OPTRC=[];

count=1;

for p=0.05:0.05:0.5

    MissM=MissM_all(:,:,:,count);
    MissM_int=uint8(MissM*255);

    Mask = uint8(zeros(n_o));
    Mask(MissM_int~=0) = 1;

    Mask_int=Mask;

    MissM=double(MissM_int)/255;
    Mask=double(Mask_int);
    
    MissM3=reshape(MissM,n_o(1),n_o(2),[]);
    Mask3=reshape(Mask,n_o(1),n_o(2),[]);

    disp(['I_OPTRC_old PSNR: ' num2str(psnr(I_OPTRC_all_old(:,:,:,count),I)) ' time:' num2str(time_OPTRC_old(count))]);
    
    %% OPTRC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    option=get_option_default(MissM_int,Mask_int,'OPTRC',opt);
    option.buffersize = 1;
    t1=tic;
    I_OPTRC=OPTRC(MissM_int,Mask_int,[],option);
    I_OPTRC=(1-Mask).*I_OPTRC+Mask.*MissM;
    time_OPTRC=[time_OPTRC toc(t1)];
    disp(['I_OPTRC PSNR: ' num2str(psnr(I_OPTRC,I)) ' time:' num2str(time_OPTRC(end))]);
  
    I_OPTRC_all=cat(4,I_OPTRC_all,I_OPTRC);
    
    count=count+1;
    
end

save results\image_p_new.mat









