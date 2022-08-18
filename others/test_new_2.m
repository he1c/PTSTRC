clear

addpath(genpath(pwd))

opt.framenum=1;

imgPath = 'data\BSD\';        % 图像库路径
imgDir  = dir([imgPath '*.jpg']); % 遍历所有jpg格式文件

I_all=[];
MissM_all=[];
I_LPTR_all=[];
I_PTRC_all=[];
I_TRNN_all=[];
I_TNN_all=[];
I_TCTF_all=[];
I_TRLRF_all=[];
I_TMAC_all=[];
I_NLPTNN_all=[];
I_NLPTT_all=[];


for i = 1:1%length(imgDir)          % 遍历结构体就可以一一处理图片了

    I = imread([imgPath imgDir(i).name]); %读取每张图片
        
    disp(['filename:' imgDir(i).name]);

    if size(I,2)<size(I,1)
        I = permute(I,[2 1 3]);
    end

    I_int = imresize(I,[320 480]);
    
    I=double(I_int)/255;

    n_o=size(I_int);

    p=0.4;

    omega=rand(prod(n_o),1)<p;
    
    %G=noisemix(prod(n_o),1,0,0.000,0.25,'gaussian');
    %G=reshape(G,n_o);
    %I_int=uint8(double(I_int)+G*255);

    Mask = uint8(zeros(n_o));
    Mask(omega) = 1;

    MissM_int=Mask.*I_int;
    Mask_int=Mask;

    MissM=double(MissM_int)/255;
    Mask=double(Mask_int);

    option.I=I;
    option.r=[20,20,20];
    Ihat=TRSGD_scaled(MissM,Mask,option);
    
    %% PTRC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    opt.n_t=[4 4 4 5 4 4 5 6 3];  %320*480
    opt.eta=0.1;
    Mask_t=reshape(Mask,opt.n_t);
    MissM_t=reshape(MissM,opt.n_t);
    option=get_option_default(MissM_t,Mask_t,'PTRC',opt);
    I_hat=PTRC(MissM_t,Mask_t,[],option);
    I_PTRC=reshape(I_hat,size(I));
    disp(['PTRC SSIM: ' num2str(ssim_video(I_PTRC,I)) ]);

    %% TRNN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    option=get_option_default(MissM_t,Mask_t,'TRNN',opt);
    I_hat = TRNN(MissM_t,Mask_t,I,option);
    I_TRNN = reshape(I_hat,size(I));
    disp(['I_TRNN SSIM: ' num2str(ssim_video(I_TRNN,I)) ]);

    %% TRLRF %%%%%%%%%%%%
    option=get_option_default(MissM,Mask,'TRLRF',opt);
    [I_hat,~,~]=TRLRF(MissM,Mask,option);
    I_TRLRF=reshape(I_hat,size(I));
    disp(['I_TRLRF SSIM: ' num2str(ssim_video(I_TRLRF,I)) ]);
    
    para_TR.r=30*ones(3,1);para_TR.max_tot= 1e-5; para_TR.max_iter=100;para_TR.disp =  1;
    tic
    Utr_TR = Completion_TR(MissM, Mask, para_TR);
    time_TR=toc;
    I_hat=fullTR(Utr_TR);

   
end










