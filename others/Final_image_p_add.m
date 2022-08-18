clear

addpath(genpath(pwd))

load results\image_p_new.mat

for p=0.55:0.05:0.7

    omega=rand(prod(n_o),1)<p;

    Mask = uint8(zeros(n_o));
    Mask(omega) = 1;

    MissM_int=Mask.*I_int;
    Mask_int=Mask;

    MissM=double(MissM_int)/255;
    Mask=double(Mask_int);
    
    MissM3=reshape(MissM,n_o(1),n_o(2),[]);
    Mask3=reshape(Mask,n_o(1),n_o(2),[]);

    %% OPTRC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    option=get_option_default(MissM_int,Mask_int,'OPTRC',opt);
    option.buffersize = 1;
    t1=tic;
    I_OPTRC=OPTRC(MissM_int,Mask_int,[],option);
    I_OPTRC=(1-Mask).*I_OPTRC+Mask.*MissM;
    time_OPTRC=[time_OPTRC toc(t1)];
    disp(['I_OPTRC PSNR: ' num2str(psnr(I_OPTRC,I)) ' time:' num2str(time_OPTRC(end))]);
    
    %% OPTRC-TTNN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    option=get_option_default(MissM_int,Mask_int,'OPTRC',opt);
    option.buffersize = 1;
    option.beta=0.03;
    option.maxitr=200;
    option.stopc=1e-3;
    t1=tic;
    I_OPTRC_TTNN=OPTRC_TTNN(MissM_int,Mask_int,[],option);
    time_OPTRC_TTNN=[time_OPTRC_TTNN toc(t1)];
    disp(['I_OPTRC_TTNN SSIM: ' num2str(psnr(I_OPTRC_TTNN,I)) ' time:' num2str(time_OPTRC_TTNN(end))]);
    
    %% OPTRC-TRNN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    option=get_option_default(MissM_int,Mask_int,'OPTRC',opt);
    option.buffersize = 1;
    option.beta = 1e-5;
    option.alpha = 1e-4;
    t1=tic;
    I_OPTRC_TRNN=OPTRC_TRNN(MissM_int,Mask_int,[],option);
    time_OPTRC_TRNN=[time_OPTRC_TRNN toc(t1)];
    disp(['I_OPTRC_TRNN SSIM: ' num2str(psnr(I_OPTRC_TRNN,I)) ' time:' num2str(time_OPTRC_TRNN(end))]);

    %% PTRC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    opt.n_t=[4 4 4 5 4 4 5 6 3];  %320*480
    opt.eta=0.1;
    Mask_t=reshape(Mask,opt.n_t);
    MissM_t=reshape(MissM,opt.n_t);
    option=get_option_default(MissM_t,Mask_t,'PTRC',opt);
    option.r=(floor(p*100)+40)*ones(1,11);
    t1=tic;
    I_hat=PTRC_RW(MissM_t,Mask_t,[],option);
    I_PTRC=reshape(I_hat,size(I));
    I_PTRC=(1-Mask).*I_PTRC+Mask.*MissM;
    time_PTRC=[time_PTRC toc(t1)];
    disp(['I_PTRC SSIM: ' num2str(psnr(I_PTRC,I)) ' time:' num2str(time_PTRC(end))]);

    %% TRNN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    option=get_option_default(MissM_t,Mask_t,'TRNN',opt);
    t1=tic;
    I_hat = TRNN(MissM_t,Mask_t,I,option);
    I_TRNN = reshape(I_hat,size(I));
    time_TRNN=[time_TRNN toc(t1)];
    disp(['I_TRNN SSIM: ' num2str(psnr(I_TRNN,I)) ' time:' num2str(time_TRNN(end))]);

    %% TNN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    option=get_option_default(MissM3,Mask3,'TNN',opt);
    t1=tic;
    A = diag(sparse(double(Mask3(:)))); 
    b = MissM3(:);
    I_hat = TNN(A,b,option);
    I_TNN = reshape(I_hat,size(I));
    time_TNN=[time_TNN toc(t1)];
    disp(['I_TNN SSIM: ' num2str(psnr(I_TNN,I)) ' time:' num2str(time_TNN(end))]);

    %% TCTF %%%%%%%%%%
    option=get_option_default(MissM3,Mask3,'TCTF',opt);
    option.rank=floor(p*100);
    t1=tic;
    I_hat = TCTF(MissM3,Mask3,I,option);
    I_TCTF = reshape(I_hat,size(I));
    I_TCTF=(1-Mask).*I_TCTF+Mask.*MissM;
    time_TCTF=[time_TCTF toc(t1)];
    disp(['I_TCTF SSIM: ' num2str(psnr(I_TCTF,I)) ' time:' num2str(time_TCTF(end))]);

    %% TMAC %%%%%%%%%%%%
    option=get_option_default(MissM3,Mask3,'TMAC',opt);
    t1=tic;
    I_hat= TMAC(MissM3,Mask3,[],option);
    I_TMAC=reshape(I_hat,size(I));
    I_TMAC=(1-Mask).*I_TMAC+Mask.*MissM;
    time_TMAC=[time_TMAC toc(t1)];
    disp(['I_TMAC SSIM: ' num2str(psnr(I_TMAC,I)) ' time:' num2str(time_TMAC(end))]);

    %% NLPTNN %%%%%%%%%%%%
    option=get_option_default(MissM,Mask,'NLPTNN',opt);
    t1=tic;
    option.THRESHOLD = 29;
    I_hat=NLPTNN(MissM,Mask,[],option);
    I_NLPTNN=reshape(I_hat,size(I));
    I_NLPTNN=(1-Mask).*I_NLPTNN+Mask.*MissM;
    time_NLPTNN=[time_NLPTNN toc(t1)];
    disp(['I_NLPTNN SSIM: ' num2str(psnr(I_NLPTNN,I)) ' time:' num2str(time_NLPTNN(end))]);

    %% NLPTT %%%%%%%%%%%%
    option=get_option_default(MissM,Mask,'NLPTT',opt);
    option.K=30;
    t1=tic;
    I_NLPTT=NLPTT(MissM,Mask,I,option);
    I_NLPTT=(1-Mask).*I_NLPTT+Mask.*MissM;
    time_NLPTT=[time_NLPTT toc(t1)];
    disp(['I_NLPTT SSIM: ' num2str(psnr(I_NLPTT,I)) ' time:' num2str(time_NLPTT(end))]);
    
    %% store data
    MissM_all=cat(4,MissM_all,MissM);
    I_OPTRC_all=cat(4,I_OPTRC_all,I_OPTRC);
    I_OPTRC_TTNN_all=cat(4,I_OPTRC_TTNN_all,I_OPTRC_TTNN);
    I_OPTRC_TRNN_all=cat(4,I_OPTRC_TRNN_all,I_OPTRC_TRNN);
    I_PTRC_all=cat(4,I_PTRC_all,I_PTRC);
    I_TRNN_all=cat(4,I_TRNN_all,I_TRNN);
    I_TNN_all=cat(4,I_TNN_all,I_TNN);
    I_TCTF_all=cat(4,I_TCTF_all,I_TCTF);
    I_TMAC_all=cat(4,I_TMAC_all,I_TMAC);
    I_NLPTNN_all=cat(4,I_NLPTNN_all,I_NLPTNN);
    I_NLPTT_all=cat(4,I_NLPTT_all,I_NLPTT);
    
end

save results\image_p_add.mat









