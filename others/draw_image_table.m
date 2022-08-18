PSNR_all=[];
PSNR_SSIM_all=[];

ind=[7,23,27,28];
ind2=[0.05,0.1,0.2];

for kkk=4:-1:1
    
    for cc=1:1:3
        
        load(['bsd_' num2str(ind2(cc)) '.mat'])
        
        if kkk>0
            load(['bsd_' num2str(ind2(cc)) '_LPTR_all2.mat'])
        else
            load(['bsd_' num2str(ind2(cc)) '_LPTR_178054.mat'])
        end
        
        PSNR=[];
        indi=ind(kkk);
        PSNR=[PSNR psnr(I_TMAC_all(:,:,:,indi),I_all(:,:,:,indi)) psnr(I_TNN_all(:,:,:,indi),I_all(:,:,:,indi)) psnr(I_TCTF_all(:,:,:,indi),I_all(:,:,:,indi)) psnr(I_TRNN_all(:,:,:,indi),I_all(:,:,:,indi)) psnr(I_PTRC_all(:,:,:,indi),I_all(:,:,:,indi))];
        PSNR=[PSNR psnr(I_NLPTNN_all(:,:,:,indi),I_all(:,:,:,indi)) psnr(I_NLPTT_all(:,:,:,indi),I_all(:,:,:,indi)) psnr(I_LPTR_NN_all(:,:,:,kkk),I_all(:,:,:,indi)) psnr(I_LPTR_MF_all(:,:,:,kkk),I_all(:,:,:,indi)) psnr(I_LPTR_ND_all(:,:,:,kkk),I_all(:,:,:,indi)) psnr(I_LPTR_all(:,:,:,kkk),I_all(:,:,:,indi))];
        SSIM=[];
        SSIM=[SSIM ssim(I_TMAC_all(:,:,:,indi),I_all(:,:,:,indi)) ssim(I_TNN_all(:,:,:,indi),I_all(:,:,:,indi)) ssim(I_TCTF_all(:,:,:,indi),I_all(:,:,:,indi)) ssim(I_TRNN_all(:,:,:,indi),I_all(:,:,:,indi)) ssim(I_PTRC_all(:,:,:,indi),I_all(:,:,:,indi))];
        SSIM=[SSIM ssim(I_NLPTNN_all(:,:,:,indi),I_all(:,:,:,indi)) ssim(I_NLPTT_all(:,:,:,indi),I_all(:,:,:,indi)) ssim(I_LPTR_NN_all(:,:,:,kkk),I_all(:,:,:,indi)) ssim(I_LPTR_MF_all(:,:,:,kkk),I_all(:,:,:,indi)) ssim(I_LPTR_ND_all(:,:,:,kkk),I_all(:,:,:,indi)) ssim(I_LPTR_all(:,:,:,kkk),I_all(:,:,:,indi))];

        PSNR_SSIM_all=[PSNR_SSIM_all;roundn(PSNR,-2);roundn(SSIM,-4)];
    end
end





    