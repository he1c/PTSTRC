clear

addpath(genpath(pwd))

%% start parallel computation
parfor i=1
    aaaa=1;
end

%% observation rate
p=0.2;

%% read video
video_num=1;
v = VideoReader(['data\Moments_in_Time_Raw\D' num2str(video_num) '.mp4']);
I_int=uint8(zeros(v.Height,v.Width,3,v.NumFrames));
count=1;
for i=1:v.NumFrames
    I_int(:,:,:,count)=imresize(read(v,i),1);
    count=count+1;
end

if size(I_int,1)>256
    I_int=imresize(I_int,256/size(I_int,1));
end

frame_num=min(50,size(I_int,4));
opt.framenum=frame_num;

I_int=I_int(:,:,:,1:frame_num);

n_o=size(I_int);

fprintf('number: %d, height: %d, width %d\n',nn,n_o(1),n_o(2));

omega=rand(prod(n_o),1)<p;

Mask = uint8(zeros(n_o));
Mask(omega) = 1;

MissM_int=Mask.*I_int;
Mask_int=Mask;

MissM=double(MissM_int)/255;
Mask=double(Mask_int);
I=double(I_int)/255;

MissM3=reshape(MissM,n_o(1),n_o(2),[]);
Mask3=reshape(Mask,n_o(1),n_o(2),[]);

%% PTSTRC %%%%%%%%%%
option=get_option_default(MissM_int,Mask_int,'OPTRC',opt);
t1=tic;
I_OPTRC=OPTRC_new(MissM_int,Mask_int,I,option);
I_OPTRC=(1-Mask).*I_OPTRC+Mask.*MissM;
time(1)=toc(t1);
disp(['I_OPTRC PSNR: ' num2str(psnr(I_OPTRC,I)) ' time:' num2str(time(1))]);

%% PTRC %%%%%%%%%%
opt.n_t=n_o;
opt.eta=0.2;
option=get_option_default(MissM,Mask,'PTRC',opt);
t1=tic;
I_hat=PTRC_RW(MissM,Mask,[],option);
I_PTRC=reshape(I_hat,size(I));
I_PTRC=(1-Mask).*I_PTRC+Mask.*MissM;
time(2)=toc(t1);
disp(['I_PTRC PSNR: ' num2str(psnr(I_PTRC,I)) ' time:' num2str(time(2))]);

%% TRNN %%%%%%%%%%
option=get_option_default(MissM,Mask,'TRNN',opt);
t1=tic;
I_hat = TRNN(MissM,Mask,I,option);
I_TRNN = reshape(I_hat,size(I));
time(3)=toc(t1);
disp(['I_TRNN PSNR: ' num2str(psnr(I_TRNN,I)) ' time:' num2str(time(3))]);

%% TCTF %%%%%%%%%%
I_TCTF=imresize(I_TCTF,[n_o(1) n_o(2)]);
time(4)=time(4)*0.6;
disp(['I_TCTF PSNR: ' num2str(psnr(I_TCTF,I)) ' time:' num2str(time(4))]);

%% TMAC %%%%%%%%%%%%
I_TMAC=imresize(I_TMAC,[n_o(1) n_o(2)]);
time(5)=time(5)*0.6;
disp(['I_TMAC PSNR: ' num2str(psnr(I_TMAC,I)) ' time:' num2str(time(5))]);

%% OLSTEC %%%%%%%%%%%%
I_OLSTEC=imresize(I_OLSTEC,[n_o(1) n_o(2)]);
time(6)=time(6)*0.6;
disp(['I_OLSTEC PSNR: ' num2str(psnr(I_OLSTEC,I)) ' time:' num2str(time(6))]);

%% TNN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
option=get_option_default(MissM3,Mask3,'TNN',opt);
t1=tic;
A = diag(sparse(double(Mask3(:)))); 
b = MissM3(:);
I_hat = TNN(A,b,option);
I_TNN = reshape(I_hat,size(I));
time(7)=toc(t1);
disp(['I_TNN PSNR: ' num2str(psnr(I_TNN,I)) ' time:' num2str(time(7))]);

%% TOUCAN %%%%%%%%%
t1=tic;
I_hat = TOUCAN(MissM3,Mask3,[],[]);
I_TOUCAN=reshape(I_hat,size(I));
I_TOUCAN=(1-Mask).*I_TOUCAN+Mask.*MissM;
time(8)=toc(t1);
disp(['I_TOUCAN PSNR: ' num2str(psnr(I_TOUCAN,I)) ' time:' num2str(time(8))]);

%% NLPTNN %%%%%%%%%%%%
option=get_option_default(MissM,Mask,'NLPTNN',opt);
t1=tic;
I_hat=NLPTNN(MissM,Mask,[],option);
I_NLPTNN=reshape(I_hat,size(I));
I_NLPTNN=(1-Mask).*I_NLPTNN+Mask.*MissM;
time(9)=toc(t1);
disp(['I_NLPTNN PSNR: ' num2str(psnr(I_NLPTNN,I)) ' time:' num2str(time(9))]);

%% NLPTT %%%%%%%%%%%%
option=get_option_default(MissM,Mask,'NLPTT',opt);
option.debug=0;
t1=tic;
I_NLPTT=NLPTT(MissM,Mask,I,option);
I_NLPTT=(1-Mask).*I_NLPTT+Mask.*MissM;
time(10)=toc(t1);
disp(['I_NLPTT PSNR: ' num2str(psnr(I_NLPTT,I)) ' time:' num2str(time(10))]);











