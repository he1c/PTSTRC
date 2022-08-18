clear

addpath(genpath(pwd))

load data\stefan.mat

frame_num=1;
opt.framenum=frame_num;

I_int=I(:,:,:,1:frame_num);

n_o=size(I_int);

p=0.3;

omega=rand(prod(n_o),1)<p;
%G=noisemix(prod(n_o),1,0,0.0,0.25,'gaussian');
%G=reshape(G,n_o);
%I_int=I_int+uint8(G*255);

Mask = uint8(zeros(n_o));
Mask(omega) = 1;

MissM_int=Mask.*I_int;
Mask_int=Mask;

MissM=double(MissM_int)/255;
Mask=double(Mask_int);
I=double(I_int)/255;

%% LPTR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
option=get_option_default(MissM_int,Mask_int,'LPTR',opt);
option.adjnum=0;
option.TRrank = 4;
option.thre_TRNN=100;
I_LPTR_MF=LPTR_mix_test(MissM_int,Mask_int,I_int,option);
   











