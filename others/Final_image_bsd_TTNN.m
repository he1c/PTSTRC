clear

addpath(genpath(pwd))

opt.framenum=1;

imgPath = 'data\BSD\';        % 图像库路径
imgDir  = dir([imgPath '*.jpg']); % 遍历所有jpg格式文件

I_LPTR_all=[];
I_LPTR_MF_all=[];
I_LPTR_NN_all=[];
I_LPTR_ND_all=[];

ind=[7,23,27,28];

% for kk = 1:1       % 遍历结构体就可以一一处理图片了
%     
%     i=ind(kk);
% 
%     I = imread([imgPath imgDir(i).name]); %读取每张图片
%         
%     disp(['filename:' imgDir(i).name]);
% 
%     if size(I,2)<size(I,1)
%         I = permute(I,[2 1 3]);
%     end
% 
%     I_int = imresize(I,[320 480]);
% 
%     n_o=size(I_int);
% 
%     p=0.05;
% 
%     omega=rand(prod(n_o),1)<p;
% 
%     Mask = uint8(zeros(n_o));
%     Mask(omega) = 1;
% 
%     MissM_int=Mask.*I_int;
%     Mask_int=Mask;
% 
%     MissM=double(MissM_int)/255;
%     Mask=double(Mask_int);
%     I=double(I_int)/255;
% 
%     MissM3=reshape(MissM,n_o(1),n_o(2),[]);
%     Mask3=reshape(Mask,n_o(1),n_o(2),[]);
% 
%     %% LPTR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     option=get_option_default(MissM_int,Mask_int,'LPTR',opt);
%     option.adjnum=0;
%     option.TRrank = 4;
%     if p>0.2
%         option.thre_TRNN=10000000;
%     end
%     I_LPTR=LPTR_mix(MissM_int,Mask_int,[],option);
% 
%     I_LPTR_all=cat(4,I_LPTR_all,I_LPTR);
%     
%     %% LPTR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     option=get_option_default(MissM_int,Mask_int,'LPTR',opt);
%     option.adjnum=0;
%     option.TRrank = 4;
%     option.thre_TRNN=10000000;
%     I_LPTR_NN=LPTR_mix(MissM_int,Mask_int,[],option);
% 
%     I_LPTR_NN_all=cat(4,I_LPTR_NN_all,I_LPTR_NN);
%     
%     %% LPTR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     option=get_option_default(MissM_int,Mask_int,'LPTR',opt);
%     option.adjnum=0;
%     option.TRrank = 4;
%     option.thre_TRNN=-1;
%     I_LPTR_MF=LPTR_mix(MissM_int,Mask_int,[],option);
% 
%     I_LPTR_MF_all=cat(4,I_LPTR_MF_all,I_LPTR_MF);
%    
%     %% LPTR-no dilation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     option=get_option_default(MissM_int,Mask_int,'LPTR',opt);
%     option.adjnum=0;
%     option.TRrank = 4;
%     option.se=strel('square',1);
%     I_LPTR_ND=LPTR_mix(MissM_int,Mask_int,[],option);
% 
%     I_LPTR_ND_all=cat(4,I_LPTR_ND_all,I_LPTR_ND);
%     
% end
% 
% save results\bsd_0.05_LPTR_178054.mat
% 
% 
% 
% I_LPTR_all=[];
% I_LPTR_MF_all=[];
% I_LPTR_NN_all=[];
% I_LPTR_ND_all=[];
% 
% 
% for kk = 1:1       % 遍历结构体就可以一一处理图片了
%     
%     i=ind(kk);
% 
%     I = imread([imgPath imgDir(i).name]); %读取每张图片
%         
%     disp(['filename:' imgDir(i).name]);
% 
%     if size(I,2)<size(I,1)
%         I = permute(I,[2 1 3]);
%     end
% 
%     I_int = imresize(I,[320 480]);
% 
%     n_o=size(I_int);
% 
%     p=0.1;
% 
%     omega=rand(prod(n_o),1)<p;
% 
%     Mask = uint8(zeros(n_o));
%     Mask(omega) = 1;
% 
%     MissM_int=Mask.*I_int;
%     Mask_int=Mask;
% 
%     MissM=double(MissM_int)/255;
%     Mask=double(Mask_int);
%     I=double(I_int)/255;
% 
%     MissM3=reshape(MissM,n_o(1),n_o(2),[]);
%     Mask3=reshape(Mask,n_o(1),n_o(2),[]);
% 
%     %% LPTR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     option=get_option_default(MissM_int,Mask_int,'LPTR',opt);
%     option.adjnum=0;
%     option.TRrank = 5;
%     if p>0.2
%         option.thre_TRNN=10000000;
%     end
%     I_LPTR=LPTR_mix(MissM_int,Mask_int,[],option);
% 
%     I_LPTR_all=cat(4,I_LPTR_all,I_LPTR);
%     
%     %% LPTR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     option=get_option_default(MissM_int,Mask_int,'LPTR',opt);
%     option.adjnum=0;
%     option.TRrank = 5;
%     option.thre_TRNN=10000000;
%     I_LPTR_NN=LPTR_mix(MissM_int,Mask_int,[],option);
% 
%     I_LPTR_NN_all=cat(4,I_LPTR_NN_all,I_LPTR_NN);
%     
%     %% LPTR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     option=get_option_default(MissM_int,Mask_int,'LPTR',opt);
%     option.adjnum=0;
%     option.TRrank = 5;
%     option.thre_TRNN=-1;
%     I_LPTR_MF=LPTR_mix(MissM_int,Mask_int,[],option);
% 
%     I_LPTR_MF_all=cat(4,I_LPTR_MF_all,I_LPTR_MF);
%    
%     %% LPTR-no dilation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     option=get_option_default(MissM_int,Mask_int,'LPTR',opt);
%     option.adjnum=0;
%     option.TRrank = 5;
%     option.se=strel('square',1);
%     if p>0.2
%         option.thre_TRNN=10000000;
%     end
%     I_LPTR_ND=LPTR_mix(MissM_int,Mask_int,[],option);
% 
%     I_LPTR_ND_all=cat(4,I_LPTR_ND_all,I_LPTR_ND);
%     
% end
% 
% save results\bsd_0.1_LPTR_178054.mat
% 


I_LPTR_all=[];
I_LPTR_MF_all=[];
I_LPTR_NN_all=[];
I_LPTR_ND_all=[];

for kk = 1:4       % 遍历结构体就可以一一处理图片了
    
    i=ind(kk);

    I = imread([imgPath imgDir(i).name]); %读取每张图片
        
    disp(['filename:' imgDir(i).name]);

    if size(I,2)<size(I,1)
        I = permute(I,[2 1 3]);
    end

    I_int = imresize(I,[320 480]);

    n_o=size(I_int);

    p=0.2;

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
    option.TRrank = 6;
    I_LPTR=LPTR_mix(MissM_int,Mask_int,[],option);

    I_LPTR_all=cat(4,I_LPTR_all,I_LPTR);
    
    %% LPTR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    option=get_option_default(MissM_int,Mask_int,'LPTR',opt);
    option.adjnum=0;
    option.TRrank = 6;
    option.thre_TRNN=10000000;
    I_LPTR_NN=LPTR_mix(MissM_int,Mask_int,[],option);

    I_LPTR_NN_all=cat(4,I_LPTR_NN_all,I_LPTR_NN);
    
    %% LPTR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    option=get_option_default(MissM_int,Mask_int,'LPTR',opt);
    option.adjnum=0;
    option.TRrank = 6;
    option.thre_TRNN=-1;
    I_LPTR_MF=LPTR_mix(MissM_int,Mask_int,[],option);

    I_LPTR_MF_all=cat(4,I_LPTR_MF_all,I_LPTR_MF);
   
    %% LPTR-no dilation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    option=get_option_default(MissM_int,Mask_int,'LPTR',opt);
    option.adjnum=0;
    option.TRrank = 6;
    option.se=strel('square',1);
    I_LPTR_ND=LPTR_mix(MissM_int,Mask_int,[],option);

    I_LPTR_ND_all=cat(4,I_LPTR_ND_all,I_LPTR_ND);
    
end

save results\bsd_0.2_LPTR_all.mat









