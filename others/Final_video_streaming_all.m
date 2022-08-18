clear

addpath(genpath(pwd))

load data\stefan.mat

frame_num=30;
opt.framenum=frame_num;

I_int=I(:,:,:,1:frame_num);

% I_int=imresize(I_int,0.6);

n_o=size(I_int);

% %%%%%%%%%%%%
% Ig=zeros(n_o(1),n_o(2),n_o(4));
% for i=1:1:n_o(4)
%     Ig(:,:,i)=rgb2gray(I(:,:,:,i));
% end
% I_int=uint8(Ig);
% n_o=size(I_int);
% %%%%%%%%%%%%

p=0.5;

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

time=[];

% %% LPTR online %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% option=get_option_default(MissM_int,Mask_int,'OPTRC',opt);
% option.TRrank = 4;
% tic
% I_LPTR_online=OPTRC(MissM_int,Mask_int,[],option);
% time=[time toc];

%% TOUCAN %%%%%%%%%
tic
I_hat = TOUCAN(MissM3,Mask3,[],[]);
I_TOUCAN=reshape(I_hat,size(I));
time=[time toc];

%% OLSTEC %%%%%%%%%%%%
option=get_option_default(MissM3,Mask,'OLSTEC',opt);
tic
[~,~, sub_infos_olstec] = OLSTEC(MissM3, Mask3, option);
I_OLSTEC=reshape(sub_infos_olstec.L,size(I));
time=[time toc];

%save results\nba_stream.mat









