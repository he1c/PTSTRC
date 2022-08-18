clear

addpath(genpath(pwd))

for p=0.4

    load(['results\image_foreman_randomstripe_' num2str(p) '.mat']);

%     parfor i=1:8
%         aa=1;
%     end
%     
%     opt.framenum=1;
%     MissM_int=MissM_int(:,:,:,25);
%     Mask_int=Mask_int(:,:,:,25);
%     I_int=I_int(:,:,:,25);
% 
%     MissM=double(MissM_int)/255;
%     Mask=double(Mask_int);
%     I=double(I_int)/255;
% 
%     MissM3=reshape(MissM,n_o(1),n_o(2),[]);
%     Mask3=reshape(Mask,n_o(1),n_o(2),[]);
% 
%     %% LPTR online%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     option=get_option_default(MissM_int,Mask_int,'OPTRC',opt);
%     t1=tic;
%     I_OPTRC=OPTRC(MissM_int,Mask_int,I,option);
%     I_OPTRC=(1-Mask).*I_OPTRC+Mask.*MissM;
%     time=[time toc(t1)];
%     disp(['I_OPTRC SSIM: ' num2str(psnr(I_OPTRC,I)) ' time:' num2str(time(end))]);
%     
%     %% OPTRC-TTNN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     option=get_option_default(MissM_int,Mask_int,'OPTRC',opt);
%     option.beta=0.03;
%     option.maxitr=200;
%     option.stopc=1e-3;
%     t1=tic;
%     I_OPTRC_TTNN=OPTRC_TTNN(MissM_int,Mask_int,[],option);
%     I_OPTRC_TTNN=(1-Mask).*I_OPTRC_TTNN+Mask.*MissM;
%     time=[time toc(t1)];
%     disp(['I_OPTRC_TTNN SSIM: ' num2str(psnr(I_OPTRC_TTNN,I)) ' time:' num2str(time(end))]);
% 
%     %% OPTRC-TRNN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     option=get_option_default(MissM_int,Mask_int,'OPTRC',opt);
%     option.beta = 1e-5;
%     option.alpha = 1e-4;
%     t1=tic;
%     I_OPTRC_TRNN=OPTRC_TRNN(MissM_int,Mask_int,[],option);
%     I_OPTRC_TRNN=(1-Mask).*I_OPTRC_TRNN+Mask.*MissM;
%     time=[time toc(t1)];
%     disp(['I_OPTRC_TRNN SSIM: ' num2str(psnr(I_OPTRC_TRNN,I)) ' time:' num2str(time(end))]);
% 
%     %% PTRC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     opt.n_t=[2 4 4 9 2 4 4 11 3];
%     opt.eta=0.1;
%     Mask_t=reshape(Mask,opt.n_t);
%     MissM_t=reshape(MissM,opt.n_t);
%     option=get_option_default(MissM_t,Mask_t,'PTRC',opt);
%     option.r=40*ones(1,9);
%     t1=tic;
%     I_hat=PTRC_RW(MissM_t,Mask_t,[],option);
%     I_PTRC=reshape(I_hat,size(I));
%     time=[time toc(t1)];
%     disp(['I_PTRC SSIM: ' num2str(psnr(I_PTRC,I)) ' time:' num2str(time(end))]);
% 
%     %% TRNN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     option=get_option_default(MissM_t,Mask_t,'TRNN',opt);
%     t1=tic;
%     I_hat = TRNN(MissM_t,Mask_t,I,option);
%     I_TRNN = reshape(I_hat,size(I));
%     time=[time toc(t1)];
%     disp(['I_TRNN SSIM: ' num2str(psnr(I_TRNN,I)) ' time:' num2str(time(end))]);
% 
%     %% TNN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     option=get_option_default(MissM3,Mask3,'TNN',opt);
%     t1=tic;
%     A = diag(sparse(double(Mask3(:)))); 
%     b = MissM3(:);
%     I_hat = TNN(A,b,option);
%     I_TNN = reshape(I_hat,size(I));
%     time=[time toc(t1)];
%     disp(['I_TNN SSIM: ' num2str(psnr(I_TNN,I)) ' time:' num2str(time(end))]);
% 
%     %% TCTF %%%%%%%%%%
%     option=get_option_default(MissM3,Mask3,'TCTF',opt);
%     option.rank   = 5;
%     t1=tic;
%     I_hat = TCTF(MissM3,Mask3,I,option);
%     I_TCTF = reshape(I_hat,size(I));
%     I_TCTF=(1-Mask).*I_TCTF+Mask.*MissM;
%     time=[time toc(t1)];
%     disp(['I_TCTF SSIM: ' num2str(psnr(I_TCTF,I)) ' time:' num2str(time(end))]);
% 
%     %% TMAC %%%%%%%%%%%%
%     option=get_option_default(MissM3,Mask3,'TMAC',opt);
%     option.rank_min = 5*ones(1,3);
%     option.rank_max = 5*ones(1,3);
%     t1=tic;
%     I_hat= TMAC(MissM3,Mask3,[],option);
%     I_TMAC=reshape(I_hat,size(I));
%     I_TMAC=(1-Mask).*I_TMAC+Mask.*MissM;
%     time=[time toc(t1)];
%     disp(['I_TMAC SSIM: ' num2str(psnr(I_TMAC,I)) ' time:' num2str(time(end))]);
% 
%     %% TOUCAN %%%%%%%%%
%     t1=tic;
%     I_hat = TOUCAN(MissM3,Mask3,[],[]);
%     I_TOUCAN=reshape(I_hat,size(I));
%     I_TOUCAN=(1-Mask).*I_TOUCAN+Mask.*MissM;
%     time=[time toc(t1)];
%     disp(['I_TOUCAN SSIM: ' num2str(psnr(I_TOUCAN,I)) ' time:' num2str(time(end))]);
% 
%     %% OLSTEC %%%%%%%%%%%%
%     option=get_option_default(MissM3,Mask,'OLSTEC',opt);
%     t1=tic;
%     [~,~, sub_infos_olstec] = OLSTEC(MissM3, Mask3, option);
%     I_OLSTEC=reshape(sub_infos_olstec.L,size(I));
%     I_OLSTEC=(1-Mask).*I_OLSTEC+Mask.*MissM;
%     time=[time toc(t1)];
%     disp(['I_OLSTEC SSIM: ' num2str(psnr(I_OLSTEC,I)) ' time:' num2str(time(end))]);

    %% NLPTNN %%%%%%%%%%%%
    option=get_option_default(MissM,Mask,'NLPTNN',opt);
    t1=tic;
    I_hat=NLPTNN(MissM,Mask,[],option);
    I_NLPTNN=reshape(I_hat,size(I));
    I_NLPTNN=(1-Mask).*I_NLPTNN+Mask.*MissM;
    time=[time toc(t1)];
    disp(['I_NLPTNN SSIM: ' num2str(psnr(I_NLPTNN,I)) ' time:' num2str(time(end))]);

    %% NLPTT %%%%%%%%%%%%
    option=get_option_default(MissM,Mask,'NLPTT',opt);
    t1=tic;
    I_NLPTT=NLPTT(MissM,Mask,I,option);
    I_NLPTT=(1-Mask).*I_NLPTT+Mask.*MissM;
    time=[time toc(t1)];
    disp(['I_NLPTT SSIM: ' num2str(psnr(I_NLPTT,I)) ' time:' num2str(time(end))]);

    save(['results\image_foreman_randomstripe_' num2str(p) '.mat']);

end






