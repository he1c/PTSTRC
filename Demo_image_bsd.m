clear

addpath(genpath(pwd))

opt.framenum=1;

%% start parallel computation
parfor i=1
    aaaa=1;
end

imgPath = 'data\BSD\';        % 图像库路径
imgDir  = dir([imgPath '*.jpg']); % 遍历所有jpg格式文件

I_all=[];
MissM_all=[];
I_PMTRSSD_all=[];
I_PMTTNN_all=[];
I_PMTRNN_all=[];
I_PTRC_all=[];
I_TRNN_all=[];
I_TNN_all=[];
I_TCTF_all=[];
I_TMAC_all=[];
I_NLPTNN_all=[];
I_NLPTT_all=[];

time_PMTRSSD=[];
time_PMTTNN=[];
time_PMTRNN=[];
time_PTRC=[];
time_TRNN=[];
time_TNN=[];
time_TCTF=[];
time_TMAC=[];
time_NLPTNN=[];
time_NLPTT=[];


for i = 1:length(imgDir)          % 遍历结构体就可以一一处理图片了

    I = imread([imgPath imgDir(i).name]); %读取每张图片
        
    disp(['filename:' imgDir(i).name]);

    if size(I,2)<size(I,1)
        I = permute(I,[2 1 3]);
    end

    I_int = imresize(I,[320 480]);
    
    I=double(I_int)/255;

    n_o=size(I_int);

    p=0.2;

    omega=rand(prod(n_o),1)<p;

    Mask = uint8(zeros(n_o));
    Mask(omega) = 1;

    MissM_int=Mask.*I_int;
    Mask_int=Mask;

    MissM=double(MissM_int)/255;
    Mask=double(Mask_int);
    
    MissM3=reshape(MissM,n_o(1),n_o(2),[]);
    Mask3=reshape(Mask,n_o(1),n_o(2),[]);

    %% PMTRSSD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    option=get_option_default(MissM_int,Mask_int,'OPTRC',opt);
    t1=tic;
    I_PMTRSSD=PTSTRC(MissM_int,Mask_int,[],option);
    I_PMTRSSD=(1-Mask).*I_PMTRSSD+Mask.*MissM;
    time_PMTRSSD=[time_PMTRSSD toc(t1)];
    disp(['I_OPTRC PSNR: ' num2str(psnr(I_PMTRSSD,I)) ' time:' num2str(time_PMTRSSD(end))]);
    
   %% PMTRNN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    option=get_option_default(MissM_int,Mask_int,'OPTRC',opt);
    option.buffersize = 1;
    option.beta = 1e-5;
    option.alpha = 1e-4;
    t1=tic;
    I_PMTRNN=PMTRNN(MissM_int,Mask_int,[],option);
    I_PMTRNN=(1-Mask).*I_PMTRNN+Mask.*MissM;
    time_PMTRNN=[time_PMTRNN toc(t1)];
    disp(['I_OPTRC_TRNN PSNR: ' num2str(psnr(I_PMTRNN,I)) ' time:' num2str(time_PMTRNN(end))]);
      
    %% PMTTNN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    option=get_option_default(MissM_int,Mask_int,'OPTRC',opt);
    option.buffersize = 1;
    option.beta=0.03;
    option.maxitr=200;
    option.stopc=1e-3;
    t1=tic;
    I_PMTTNN=PMTTNN(MissM_int,Mask_int,[],option);
    I_PMTTNN=(1-Mask).*I_PMTTNN+Mask.*MissM;
    time_PMTTNN=[time_PMTTNN toc(t1)];
    disp(['I_OPTRC_TTNN PSNR: ' num2str(psnr(I_PMTTNN,I)) ' time:' num2str(time_PMTTNN(end))]);

    %% PTRC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    opt.n_t=[4 4 4 5 4 4 5 6 3];  %320*480
    opt.eta=0.1;
    Mask_t=reshape(Mask,opt.n_t);
    MissM_t=reshape(MissM,opt.n_t);
    option=get_option_default(MissM_t,Mask_t,'PTRC',opt);
    option.r=60*ones(1,11);
    t1=tic;
    I_hat=PTRC_RW(MissM_t,Mask_t,[],option);
    I_PTRC=reshape(I_hat,size(I));
    time_PTRC=[time_PTRC toc(t1)];
    disp(['I_PTRC PSNR: ' num2str(psnr(I_PTRC,I)) ' time:' num2str(time_PTRC(end))]);

    %% TRNN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    option=get_option_default(MissM_t,Mask_t,'TRNN',opt);
    t1=tic;
    I_hat = TRNN(MissM_t,Mask_t,I,option);
    I_TRNN = reshape(I_hat,size(I));
    time_TRNN=[time_TRNN toc(t1)];
    disp(['I_TRNN PSNR: ' num2str(psnr(I_TRNN,I)) ' time:' num2str(time_TRNN(end))]);

    %% TNN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    option=get_option_default(MissM3,Mask3,'TNN',opt);
    t1=tic;
    A = diag(sparse(double(Mask3(:)))); 
    b = MissM3(:);
    I_hat = TNN(A,b,option);
    I_TNN = reshape(I_hat,size(I));
    time_TNN=[time_TNN toc(t1)];
    disp(['I_TNN PSNR: ' num2str(psnr(I_TNN,I)) ' time:' num2str(time_TNN(end))]);

    %% TCTF %%%%%%%%%%
    option=get_option_default(MissM3,Mask3,'TCTF',opt);
    t1=tic;
    I_hat = TCTF(MissM3,Mask3,I,option);
    I_TCTF = reshape(I_hat,size(I));
    time_TCTF=[time_TCTF toc(t1)];
    disp(['I_TCTF PSNR: ' num2str(psnr(I_TCTF,I)) ' time:' num2str(time_TCTF(end))]);

    %% TMAC %%%%%%%%%%%%
    option=get_option_default(MissM3,Mask3,'TMAC',opt);
    t1=tic;
    I_hat= TMAC(MissM3,Mask3,[],option);
    I_TMAC=reshape(I_hat,size(I));
    time_TMAC=[time_TMAC toc(t1)];
    disp(['I_TMAC PSNR: ' num2str(psnr(I_TMAC,I)) ' time:' num2str(time_TMAC(end))]);
    
%     %% TOUCAN %%%%%%%%%
%     I_hat = TOUCAN(MissM3,Mask3,[],[]);
%     I_TOUCAN=reshape(I_hat,size(I));
%     disp(['I_TOUCAN PSNR: ' num2str(psnr(I_TOUCAN,I)) ]);

    %% OLSTEC %%%%%%%%%%%%
    option=get_option_default(MissM3,Mask,'OLSTEC',opt);
    [~,~, sub_infos_olstec] = OLSTEC(MissM3, Mask3, option);
    I_OLSTEC=reshape(sub_infos_olstec.L,size(I));
    disp(['I_OLSTEC PSNR: ' num2str(psnr(I_OLSTEC,I)) ]);

    %% NLPTNN %%%%%%%%%%%%
    option=get_option_default(MissM,Mask,'NLPTNN',opt);
    t1=tic;
    I_hat=NLPTNN(MissM,Mask,[],option);
    I_NLPTNN=reshape(I_hat,size(I));
    time_NLPTNN=[time_NLPTNN toc(t1)];
    disp(['I_NLPTNN PSNR: ' num2str(psnr(I_NLPTNN,I)) ' time:' num2str(time_NLPTNN(end))]);

    %% NLPTT %%%%%%%%%%%%
    option=get_option_default(MissM_int,Mask_int,'OPTRC',opt);
    option.buffersize = 1;
    option.beta=0.1;
    option.maxitr=100;
    option.stopc=1e-3;
    option.interval=1;
    option.se=strel('square',1);
    t1=tic;
    I_NLPTT=PMTTNN(MissM_int,Mask_int,[],option);
    time_NLPTT=[time_NLPTT toc(t1)];
    disp(['I_NLPTT PSNR: ' num2str(psnr(I_NLPTT,I)) ' time:' num2str(time_NLPTT(end))]);
    
    %% store data
    I_all=cat(4,I_all,I);
    MissM_all=cat(4,MissM_all,MissM);
    I_PMTRSSD_all=cat(4,I_PMTRSSD_all,I_PMTRSSD);
    I_PMTTNN_all=cat(4,I_PMTTNN_all,I_PMTTNN);
    I_PMTRNN_all=cat(4,I_PMTRNN_all,I_PMTRNN);
    I_PTRC_all=cat(4,I_PTRC_all,I_PTRC);
    I_TRNN_all=cat(4,I_TRNN_all,I_TRNN);
    I_TNN_all=cat(4,I_TNN_all,I_TNN);
    I_TCTF_all=cat(4,I_TCTF_all,I_TCTF);
    I_TMAC_all=cat(4,I_TMAC_all,I_TMAC);
    I_NLPTNN_all=cat(4,I_NLPTNN_all,I_NLPTNN);
    I_NLPTT_all=cat(4,I_NLPTT_all,I_NLPTT);
    I_TOUCAN_all=cat(4,I_TOUCAN_all,I_TOUCAN);
    I_OLSTEC_all=cat(4,I_OLSTEC_all,I_OLSTEC);
    
end









