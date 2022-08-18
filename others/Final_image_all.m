clear

addpath(genpath(pwd))

opt.framenum=1;

I_int=imread('data\lena.bmp');

n_o=size(I_int);

p=0.1;

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

%% LPTR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
option=get_option_default(MissM_int,Mask_int,'LPTR',opt);
option.adjnum=0;
option.TRrank = 5;
I_LPTR=LPTR_mix(MissM_int,Mask_int,[],option);

%% PTRC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%opt.n_t=[2 4 4 9 2 4 4 11 3];
%opt.n_t=[4 4 4 4 4 4 4 4 3];
opt.n_t=[4 4 4 5 4 4 5 6 3];  %320*480
Mask_t=reshape(Mask,opt.n_t);
MissM_t=reshape(MissM,opt.n_t);
option=get_option_default(MissM_t,Mask_t,'PTRC',opt);
I_hat=PTRC(MissM_t,Mask_t,[],option);
I_PTRC=reshape(I_hat,size(I));

%% TRNN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
option=get_option_default(MissM_t,Mask_t,'TRNN',opt);
I_hat = TRNN(MissM_t,Mask_t,I,option);
I_TRNN = reshape(I_hat,size(I));

%% TNN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
option=get_option_default(MissM3,Mask3,'TNN',opt);
A = diag(sparse(double(Mask3(:)))); 
b = MissM3(:);
I_hat = TNN(A,b,option);
I_TNN = reshape(I_hat,size(I));

%% TCTF %%%%%%%%%%
option=get_option_default(MissM3,Mask3,'TCTF',opt);
I_hat = TCTF(MissM3,Mask3,I,option);
I_TCTF = reshape(I_hat,size(I));

%% TRLRF %%%%%%%%%%%%
option=get_option_default(MissM,Mask,'TRLRF',opt);
[I_hat,~,~]=TRLRF(MissM,Mask,option);
I_TRLRF=reshape(I_hat,size(I));

%% TMAC %%%%%%%%%%%%
option=get_option_default(MissM3,Mask3,'TMAC',opt);
I_hat= TMAC(MissM3,Mask3,[],option);
I_TMAC=reshape(I_hat,size(I));

%% NLPTNN %%%%%%%%%%%%
option=get_option_default(MissM,Mask,'NLPTNN',opt);
I_hat=NLPTNN(MissM,Mask,[],option);
I_NLPTNN=reshape(I_hat,size(I));

%% NLPTT %%%%%%%%%%%%
option=get_option_default(MissM,Mask,'NLPTT',opt);
I_NLPTT=NLPTT(MissM,Mask,I,option);








