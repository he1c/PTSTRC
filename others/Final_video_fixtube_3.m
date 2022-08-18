clear

addpath(genpath(pwd))

for p=0.1:0.1:0.3

    load(['results\video_akiyo_fixtube_' num2str(p) '.mat']);

    parfor i=1:8
        aa=1;
    end

    %% LPTR online%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    option=get_option_default(MissM_int,Mask_int,'OPTRC',opt);
    t1=tic;
    I_OPTRC=OPTRC(MissM_int,Mask_int,I,option);
    I_OPTRC=(1-Mask).*I_OPTRC+Mask.*MissM;
    time(1)=toc(t1);
    disp(['I_OPTRC SSIM: ' num2str(psnr(I_OPTRC,I)) ' time:' num2str(time(1))]);

    %% TCTF %%%%%%%%%%
    option=get_option_default(MissM3,Mask3,'TCTF',opt);
    option.rank = 5;
    t1=tic;
    I_hat = TCTF(MissM3,Mask3,I,option);
    I_TCTF = reshape(I_hat,size(I));
    I_TCTF=(1-Mask).*I_TCTF+Mask.*MissM;
    time(4)=toc(t1);
    disp(['I_TCTF SSIM: ' num2str(psnr(I_TCTF,I)) ' time:' num2str(time(4))]);

    %% TMAC %%%%%%%%%%%%
    option=get_option_default(MissM3,Mask3,'TMAC',opt);
    option.rank_min = [10,10,10];
    option.rank_max = [10,10,10];
    t1=tic;
    I_hat= TMAC(MissM3,Mask3,[],option);
    I_TMAC=reshape(I_hat,size(I));
    I_TMAC=(1-Mask).*I_TMAC+Mask.*MissM;
    time(5)=toc(t1);
    disp(['I_TMAC SSIM: ' num2str(psnr(I_TMAC,I)) ' time:' num2str(time(5))]);

    save(['results\video_akiyo_fixtube_' num2str(p) '.mat']);

end





