load results\image_bsd_0.2_new.mat

SSIM_all=[];
PSNR_all=[];
time_all=[];

% for i=1:1:30
%     SSIM=[];
%     Mask=double(MissM_all(:,:,:,i)~=0);
%     I_OPTRC_all(:,:,:,i)=(1-Mask).*I_OPTRC_all(:,:,:,i)+Mask.*MissM_all(:,:,:,i);
%     SSIM=[SSIM ssim(I_TMAC_all(:,:,:,i),I_all(:,:,:,i)) ssim(I_TNN_all(:,:,:,i),I_all(:,:,:,i)) ssim(I_TCTF_all(:,:,:,i),I_all(:,:,:,i)) ssim(I_TRNN_all(:,:,:,i),I_all(:,:,:,i)) ssim(I_PTRC_all(:,:,:,i),I_all(:,:,:,i))];
%     SSIM=[SSIM ssim(I_NLPTNN_all(:,:,:,i),I_all(:,:,:,i)) ssim(I_NLPTT_all(:,:,:,i),I_all(:,:,:,i)) ssim(I_OPTRC_TTNN_all(:,:,:,i),I_all(:,:,:,i)) ssim(I_OPTRC_TRNN_all(:,:,:,i),I_all(:,:,:,i)) ssim(I_OPTRC_all(:,:,:,i),I_all(:,:,:,i))];
%     SSIM_all=[SSIM_all;SSIM];
% end

for i=1:1:30
    PSNR=[];
    Mask=double(MissM_all(:,:,:,i)~=0);
    I_TMAC_all(:,:,:,i)=(1-Mask).*I_TMAC_all(:,:,:,i)+Mask.*MissM_all(:,:,:,i);
    I_TCTF_all(:,:,:,i)=(1-Mask).*I_TCTF_all(:,:,:,i)+Mask.*MissM_all(:,:,:,i);
    I_PTRC_all(:,:,:,i)=(1-Mask).*I_PTRC_all(:,:,:,i)+Mask.*MissM_all(:,:,:,i);
    I_NLPTNN_all(:,:,:,i)=(1-Mask).*I_NLPTNN_all(:,:,:,i)+Mask.*MissM_all(:,:,:,i);
    I_NLPTT_all(:,:,:,i)=(1-Mask).*I_NLPTT_all(:,:,:,i)+Mask.*MissM_all(:,:,:,i);
    I_OPTRC_TTNN_all(:,:,:,i)=(1-Mask).*I_OPTRC_TTNN_all(:,:,:,i)+Mask.*MissM_all(:,:,:,i);
    I_OPTRC_TRNN_all(:,:,:,i)=(1-Mask).*I_OPTRC_TRNN_all(:,:,:,i)+Mask.*MissM_all(:,:,:,i);
    I_OPTRC_all(:,:,:,i)=(1-Mask).*I_OPTRC_all(:,:,:,i)+Mask.*MissM_all(:,:,:,i);
    PSNR=[PSNR psnr(I_TMAC_all(:,:,:,i),I_all(:,:,:,i)) psnr(I_TNN_all(:,:,:,i),I_all(:,:,:,i)) psnr(I_TCTF_all(:,:,:,i),I_all(:,:,:,i)) psnr(I_TRNN_all(:,:,:,i),I_all(:,:,:,i)) psnr(I_PTRC_all(:,:,:,i),I_all(:,:,:,i))];
    PSNR=[PSNR psnr(I_NLPTNN_all(:,:,:,i),I_all(:,:,:,i)) psnr(I_NLPTT_all(:,:,:,i),I_all(:,:,:,i)) psnr(I_OPTRC_TTNN_all(:,:,:,i),I_all(:,:,:,i)) psnr(I_OPTRC_TRNN_all(:,:,:,i),I_all(:,:,:,i)) psnr(I_OPTRC_all(:,:,:,i),I_all(:,:,:,i))];
    PSNR_all=[PSNR_all;PSNR];
end

time_all=[time_TMAC;time_TNN;time_TCTF;time_TRNN;time_PTRC;time_NLPTNN;time_NLPTT;time_OPTRC_TTNN;time_OPTRC_TRNN;time_OPTRC]';


figure,
fig=bar(1:30,PSNR_all);xlabel('image number');ylabel('PSNR(dB)');
legend(fig,'TMAC','TNN','TCTF','TRNN','PTRC','NLTNN','NLTTNN','PTTTNN','PTTRNN','PTOTRC','location','northoutside','orientation','horizontal');

figure,
fig=bar(1:30,time_all);xlabel('image number');ylabel('Average running time(s)');
legend(fig,'TMAC','TNN','TCTF','TRNN','PTRC','NLTNN','NLTTNN','PTTTNN','PTTRNN','PTOTRC','location','northoutside','orientation','horizontal');