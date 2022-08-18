clear

addpath(genpath(pwd))

load data\boxing.mat
I=imresize(X,1);
clear X
frame_num=30;
opt.framenum=frame_num;

I_int=uint8(I(:,:,:,1:frame_num)*255);

n_o=size(I_int);

for p=0.3:0.1:0.3

    
    Mask2=[];
    for i=1:1:frame_num
        Mask = uint8(ones(n_o(1),n_o(2),3));
        bknum=ceil(rand(1)*50+100);
        for j=1:1:bknum
            omega=randperm(n_o(1)-20);
            omega2=randperm(n_o(2)-300);
            if rand(1)>0
                height=ceil(rand(1)*20);
                width=ceil(rand(1)*300);
                Mask(omega:omega+height,omega2:omega2+width,1) = 0;
            end
            if rand(1)>0.5
                height=ceil(rand(1)*20);
                width=ceil(rand(1)*300);
                Mask(omega:omega+height,omega2:omega2+width,2) = 0;
            end
            if rand(1)>0.5
                height=ceil(rand(1)*20);
                width=ceil(rand(1)*300);
                Mask(omega:omega+height,omega2:omega2+width,3) = 0;
            end
        end
%             for j=1:1:5
%                 omega=randperm(n_o(2)-10);
%                 Mask(:,omega:omega+10,1,k) = 0;
%                 Mask(:,omega:omega+10,2,k) = 0;
%                 Mask(:,omega:omega+10,3,k) = 0;
%             end
        Mask2=cat(4,Mask2,Mask);
    end

%     Mask2=Mask;
%     for i=1:1:frame_num-1
%         Mask2=cat(4,Mask2,Mask);
%     end
    Mask=Mask2;

    MissM_int=Mask.*I_int;
    Mask_int=Mask;

    MissM=double(MissM_int)/255;
    Mask=double(Mask_int);
    I=double(I_int)/255;

    MissM3=reshape(MissM,n_o(1),n_o(2),[]);
    Mask3=reshape(Mask,n_o(1),n_o(2),[]);

    time=[];

    parfor k=1:8
        aa=1;
    end

%     %% LPTR online%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     option=get_option_default(MissM_int,Mask_int,'OPTRC',opt);
%     option.dim=96;
%     option.search_window_size=40;     
%     option.stride=72; 
%     option.interval=4;
%     t1=tic;
%     I_OPTRC=OPTRC(MissM_int,Mask_int,I,option);
%     I_OPTRC=(1-Mask).*I_OPTRC+Mask.*MissM;
%     time=[time toc(t1)];
%     disp(['I_OPTRC SSIM: ' num2str(psnr(I_OPTRC,I)) ' time:' num2str(time(end))]);
%     
%     %% LPTR online%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     option=get_option_default(MissM_int,Mask_int,'OPTRC',opt);
%     option.dim=96;
%     option.search_window_size=40;     
%     option.stride=48; 
%     option.interval=4;
%     t1=tic;
%     I_OPTRC2=OPTRC_block(MissM_int,Mask_int,I,option);
%     I_OPTRC2=(1-Mask).*I_OPTRC2+Mask.*MissM;
%     time=[time toc(t1)];
%     disp(['I_OPTRC SSIM: ' num2str(psnr(I_OPTRC2,I)) ' time:' num2str(time(end))]);
    
    %% LPTR online%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    option=get_option_default(MissM_int,Mask_int,'OPTRC',opt);
    option.dim=96;
    option.search_window_size=40;     
    option.stride=72; 
    option.interval=4;
    t1=tic;
    I_OPTRC3=OPTRC_new(MissM_int,Mask_int,I,option);
    I_OPTRC3=(1-Mask).*I_OPTRC3+Mask.*MissM;
    time=[time toc(t1)];
    disp(['I_OPTRC SSIM: ' num2str(psnr(I_OPTRC3,I)) ' time:' num2str(time(end))]);
% 
%     %% PTRC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %opt.n_t=[2 4 4 9 2 4 4 11 3 5 10];
%     %opt.n_t=[2 4 4 9 4 4 4 8 3 5 6];
%     opt.n_t=[4 5 6 6 4 5 8 8 3 5 6];
%     opt.eta=0.1;
%     Mask_t=reshape(Mask,opt.n_t);
%     MissM_t=reshape(MissM,opt.n_t);
%     option=get_option_default(MissM_t,Mask_t,'PTRC',opt);
%     option.r=floor(0.8*70)*ones(1,11);
%     t1=tic;
%     I_hat=PTRC_RW(MissM_t,Mask_t,[],option);
%     I_PTRC=reshape(I_hat,size(I));
%     time=[time toc(t1)];
%     disp(['I_PTRC SSIM: ' num2str(psnr(I_PTRC,I)) ' time:' num2str(time(end))]);

    %% TRNN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    option=get_option_default(MissM_t,Mask_t,'TRNN',opt);
    t1=tic;
    I_hat = TRNN(MissM_t,Mask_t,I,option);
    I_TRNN = reshape(I_hat,size(I));
    time=[time toc(t1)];
    disp(['I_TRNN SSIM: ' num2str(psnr(I_TRNN,I)) ' time:' num2str(time(end))]);

    %% TNN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    option=get_option_default(MissM3,Mask3,'TNN',opt);
    t1=tic;
    A = diag(sparse(double(Mask3(:)))); 
    b = MissM3(:);
    I_hat = TNN(A,b,option);
    I_TNN = reshape(I_hat,size(I));
    time=[time toc(t1)];
    disp(['I_TNN SSIM: ' num2str(psnr(I_TNN,I)) ' time:' num2str(time(end))]);

%     %% TCTF %%%%%%%%%%
%     option=get_option_default(MissM3,Mask3,'TCTF',opt);
%     option.rank  = 10;
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
    %% TOUCAN %%%%%%%%%
    t1=tic;
    I_hat = TOUCAN(MissM3,Mask3,[],[]);
    I_TOUCAN=reshape(I_hat,size(I));
    I_TOUCAN=(1-Mask).*I_TOUCAN+Mask.*MissM;
    time=[time toc(t1)];
    disp(['I_TOUCAN SSIM: ' num2str(psnr(I_TOUCAN,I)) ' time:' num2str(time(end))]);

    %% OLSTEC %%%%%%%%%%%%
    option=get_option_default(MissM3,Mask,'OLSTEC',opt);
    t1=tic;
    [~,~, sub_infos_olstec] = OLSTEC(MissM3, Mask3, option);
    I_OLSTEC=reshape(sub_infos_olstec.L,size(I));
    I_OLSTEC=(1-Mask).*I_OLSTEC+Mask.*MissM;
    time=[time toc(t1)];
    disp(['I_OLSTEC SSIM: ' num2str(psnr(I_OLSTEC,I)) ' time:' num2str(time(end))]);
% 
%     %% NLPTNN %%%%%%%%%%%%
%     option=get_option_default(MissM,Mask,'NLPTNN',opt);
%     t1=tic;
%     I_hat=NLPTNN(MissM,Mask,[],option);
%     I_NLPTNN=reshape(I_hat,size(I));
%     I_NLPTNN=(1-Mask).*I_NLPTNN+Mask.*MissM;
%     I_NLPTNN(I_NLPTNN>1)=1;
%     I_NLPTNN(I_NLPTNN<0)=0;
%     time=[time toc(t1)];
%     disp(['I_NLPTNN SSIM: ' num2str(psnr(I_NLPTNN,I)) ' time:' num2str(time(end))]);
% 
%     %% NLPTT %%%%%%%%%%%%
%     option=get_option_default(MissM,Mask,'NLPTT',opt);
%     t1=tic;
%     I_NLPTT=NLPTT(MissM,Mask,I,option);
%     I_NLPTT=(1-Mask).*I_NLPTT+Mask.*MissM;
%     time=[time toc(t1)];
%     disp(['I_NLPTT SSIM: ' num2str(psnr(I_NLPTT,I)) ' time:' num2str(time(end))]);

    %save(['results\video_stefan_stripe_' num2str(p) '.mat']);

end











