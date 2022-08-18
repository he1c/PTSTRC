clear

addpath(genpath(pwd))

load data\bus.mat

frame_num=50;
opt.framenum=frame_num;

I_int=I(:,:,:,1:frame_num);

n_o=size(I_int);

p=0.1;

omega=rand(prod(n_o),1)<p;

Mask = uint8(zeros(n_o));
Mask(omega) = 1;

MissM_int=Mask.*I_int;
Mask_int=Mask;

clear I_int Mask


for stride=16:4:32
    

    stride

    time=[];

    %% LPTR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    option=get_option_default(MissM_int,Mask_int,'LPTR',opt);
    option.TRrank = 5;
    option.stride=option.dim-stride;
    tStart = tic;
    [I_LPTR,time_LPTR]=LPTR_mix_time(MissM_int,Mask_int,[],option);
    time=[time toc(tStart)];

    %% LPTR online %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    option=get_option_default(MissM_int,Mask_int,'LPTR_online',opt);
    option.TRrank = 5;
    option.stride=option.dim-stride;
    tStart = tic;
    [I_LPTR_online,time_LPTR_online]=LPTR_online_time(MissM_int,Mask_int,[],option);
    time=[time toc(tStart)];

    save(['results\parameter_stride_time_' num2str(stride) '.mat'],'time_LPTR','time_LPTR_online');
    

end








