clear

addpath(genpath(pwd))

opt.framenum=1;

load results\image_p_new.mat

I_TCTF_all=[];
I_PTRC_all=[];
I_TMAC_all=[];

time_TCTF=[];
time_PTRC=[];
time_TMAC=[];

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

    %% PTRC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    opt.n_t=[4 4 4 5 4 4 5 6 3];  %320*480
    opt.eta=0.1;
    Mask_t=reshape(Mask,opt.n_t);
    MissM_t=reshape(MissM,opt.n_t);
    option=get_option_default(MissM_t,Mask_t,'PTRC',opt);
    option.r=(floor(p*150)+30)*ones(1,11);
    t1=tic;
    I_hat=PTRC_RW(MissM_t,Mask_t,[],option);
    I_PTRC=reshape(I_hat,size(I));
    I_PTRC=(1-Mask).*I_PTRC+Mask.*MissM;
    time_PTRC=[time_PTRC toc(t1)];
    disp(['I_PTRC SSIM: ' num2str(psnr(I_PTRC,I)) ' time:' num2str(time_PTRC(end))]);
    
    %% TCTF %%%%%%%%%%
    option=get_option_default(MissM3,Mask3,'TCTF',opt);
    option.rank=floor(p*200)-5;
    t1=tic;
    I_hat = TCTF(MissM3,Mask3,I,option);
    I_TCTF = reshape(I_hat,size(I));
    I_TCTF=(1-Mask).*I_TCTF+Mask.*MissM;
    time_TCTF=[time_TCTF toc(t1)];
    disp(['I_TCTF SSIM: ' num2str(psnr(I_TCTF,I)) ' time:' num2str(time_TCTF(end))]);
    
    %% TMAC %%%%%%%%%%%%
    option=get_option_default(MissM3,Mask3,'TMAC',opt);
    option.rank_min = (floor(p*200)-5)*ones(1,3);
    option.rank_max = (floor(p*200)-5)*ones(1,3);
    t1=tic;
    I_hat= TMAC(MissM3,Mask3,[],option);
    I_TMAC=reshape(I_hat,size(I));
    I_TMAC=(1-Mask).*I_TMAC+Mask.*MissM;
    time_TMAC=[time_TMAC toc(t1)];
    disp(['I_TMAC SSIM: ' num2str(psnr(I_TMAC,I)) ' time:' num2str(time_TMAC(end))]);
    
    I_TCTF_all=cat(4,I_TCTF_all,I_TCTF);
    I_PTRC_all=cat(4,I_PTRC_all,I_PTRC);
    I_TMAC_all=cat(4,I_TMAC_all,I_TMAC);
    
    count=count+1;
    
end

save results\image_p_new.mat









