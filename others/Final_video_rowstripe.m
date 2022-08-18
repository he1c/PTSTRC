clear

addpath(genpath(pwd))

load data\foreman.mat

frame_num=20;
opt.framenum=frame_num;

I_int=I(:,:,:,1:frame_num);

n_o=size(I_int);

for stripe_interval=5:2:5

Mask=[];

block_size=1;

for i=1:1:frame_num
    
    Mask_k = uint8(ones(n_o(1),n_o(2),3));
    %omega=rand(n_o(1),1)<p;
    omega=mod(i,stripe_interval)+1:stripe_interval:n_o(1);
    Mask_k(omega,:,:)=0;
    Mask=cat(4,Mask,Mask_k);
    
end

MissM_int=Mask.*I_int;
Mask_int=Mask;

MissM=double(MissM_int)/255;
Mask=double(Mask_int);
I=double(I_int)/255;

MissM3=reshape(MissM,n_o(1),n_o(2),[]);
Mask3=reshape(Mask,n_o(1),n_o(2),[]);

time=[];

parfor i=1:8
    aa=1;
end

%% LPTR online%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
option=get_option_default(MissM_int,Mask_int,'OPTRC',opt);
t1=tic;
I_OPTRC=OPTRC(MissM_int,Mask_int,I,option);
I_OPTRC=(1-Mask).*I_OPTRC+Mask.*MissM;
time=[time toc(t1)];
disp(['I_OPTRC SSIM: ' num2str(ssim_video(I_OPTRC,I)) ' time:' num2str(time(end))]);

% %% PTRC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt.n_t=[3 4 4 6 2 4 4 11 3 4 5];
opt.eta=0.1;
Mask_t=reshape(Mask,opt.n_t);
MissM_t=reshape(MissM,opt.n_t);
%option=get_option_default(MissM_t,Mask_t,'PTRC',opt);
% option.r=20*ones(1,11);
% t1=tic;
% I_hat=PTRC_RW(MissM_t,Mask_t,[],option);
% I_PTRC=reshape(I_hat,size(I));
% time=[time toc(t1)];
% disp(['I_PTRC SSIM: ' num2str(ssim_video(I_PTRC,I)) ' time:' num2str(time(end))]);
% 
% %% TRNN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% option=get_option_default(MissM_t,Mask_t,'TRNN',opt);
% t1=tic;
% I_hat = TRNN(MissM_t,Mask_t,I,option);
% I_TRNN = reshape(I_hat,size(I));
% time=[time toc(t1)];
% disp(['I_TRNN SSIM: ' num2str(ssim_video(I_TRNN,I)) ' time:' num2str(time(end))]);
   
%% TNN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
option=get_option_default(MissM3,Mask3,'TNN',opt);
t1=tic;
A = diag(sparse(double(Mask3(:)))); 
b = MissM3(:);
I_hat = TNN(A,b,option);
I_TNN = reshape(I_hat,size(I));
time=[time toc(t1)];
disp(['I_TNN SSIM: ' num2str(ssim_video(I_TNN,I)) ' time:' num2str(time(end))]);

% %% TCTF %%%%%%%%%%
% option=get_option_default(MissM3,Mask3,'TCTF',opt);
% option.rank   = 20;
% t1=tic;
% I_hat = TCTF(MissM3,Mask3,I,option);
% I_TCTF = reshape(I_hat,size(I));
% I_TCTF=(1-Mask).*I_TCTF+Mask.*MissM;
% time=[time toc(t1)];
% disp(['I_TCTF SSIM: ' num2str(ssim_video(I_TCTF,I)) ' time:' num2str(time(end))]);
% 
% %% TMAC %%%%%%%%%%%%
% option=get_option_default(MissM3,Mask3,'TMAC',opt);
% t1=tic;
% I_hat= TMAC(MissM3,Mask3,[],option);
% I_TMAC=reshape(I_hat,size(I));
% I_TMAC=(1-Mask).*I_TMAC+Mask.*MissM;
% time=[time toc(t1)];
% disp(['I_TMAC SSIM: ' num2str(ssim_video(I_TMAC,I)) ' time:' num2str(time(end))]);
% 
%% TOUCAN %%%%%%%%%
t1=tic;
I_hat = TOUCAN(MissM3,Mask3,[],[]);
I_TOUCAN=reshape(I_hat,size(I));
I_TOUCAN=(1-Mask).*I_TOUCAN+Mask.*MissM;
time=[time toc(t1)];
disp(['I_TOUCAN SSIM: ' num2str(ssim_video(I_TOUCAN,I)) ' time:' num2str(time(end))]);

%% OLSTEC %%%%%%%%%%%%
option=get_option_default(MissM3,Mask,'OLSTEC',opt);
t1=tic;
[~,~, sub_infos_olstec] = OLSTEC(MissM3, Mask3, option);
I_OLSTEC=reshape(sub_infos_olstec.L,size(I));
I_OLSTEC=(1-Mask).*I_OLSTEC+Mask.*MissM;
time=[time toc(t1)];
disp(['I_OLSTEC SSIM: ' num2str(ssim_video(I_OLSTEC,I)) ' time:' num2str(time(end))]);
% 
% %% NLPTNN %%%%%%%%%%%%
% option=get_option_default(MissM,Mask,'NLPTNN',opt);
% t1=tic;
% I_hat=NLPTNN(MissM,Mask,[],option);
% I_NLPTNN=reshape(I_hat,size(I));
% time=[time toc(t1)];
% disp(['I_NLPTNN SSIM: ' num2str(ssim_video(I_NLPTNN,I)) ' time:' num2str(time(end))]);
% 
% %% NLPTT %%%%%%%%%%%%
% option=get_option_default(MissM,Mask,'NLPTT',opt);
% option.stopc=1e-3;
% t1=tic;
% I_NLPTT=NLPTT_video(MissM,Mask,I,option);
% time=[time toc(t1)];
% disp(['I_NLPTT SSIM: ' num2str(ssim_video(I_NLPTT,I)) ' time:' num2str(time(end))]);

%save(['results\video_akiyo_diagonalstripe_' num2str(stripe_interval) '.mat']);

end






