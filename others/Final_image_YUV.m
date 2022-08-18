clear

load('results\video_flower_word.mat')
disp(['I_OPTRC PSNR: ' num2str(psnr(I_OPTRC,I)) ' time:' num2str(time(1))]);
option=get_option_default(MissM_int,Mask_int,'OPTRC',opt);
t1=tic;
I_OPTRC=OPTRC(MissM_int,Mask_int,I,option);
I_OPTRC=(1-Mask).*I_OPTRC+Mask.*MissM;
time(1)=toc(t1);
disp(['I_OPTRC2 PSNR: ' num2str(psnr(I_OPTRC,I)) ' time:' num2str(time(1))]);

save('results\video_flower_word.mat')

% for kk=0.1:0.1:0.3
% 
% load(['results\video_bus_random_' num2str(kk) '.mat'])
% disp(['I_OPTRC PSNR: ' num2str(psnr(I_OPTRC,I)) ' time:' num2str(time(1))]);
% option=get_option_default(MissM_int,Mask_int,'OPTRC',opt);
% t1=tic;
% I_OPTRC=OPTRC(MissM_int,Mask_int,I,option);
% I_OPTRC=(1-Mask).*I_OPTRC+Mask.*MissM;
% time(1)=toc(t1);
% disp(['I_OPTRC2 PSNR: ' num2str(psnr(I_OPTRC,I)) ' time:' num2str(time(1))]);
% 
% save(['results\video_bus_random_' num2str(kk) '.mat'])
% 
% end
% 
% for kk=0.1:0.1:0.3
% 
% load(['results\video_akiyo_fixtube_' num2str(kk) '.mat'])
% disp(['I_OPTRC PSNR: ' num2str(psnr(I_OPTRC,I)) ' time:' num2str(time(1))]);
% option=get_option_default(MissM_int,Mask_int,'OPTRC',opt);
% t1=tic;
% I_OPTRC=OPTRC(MissM_int,Mask_int,I,option);
% I_OPTRC=(1-Mask).*I_OPTRC+Mask.*MissM;
% time(1)=toc(t1);
% disp(['I_OPTRC PSNR: ' num2str(psnr(I_OPTRC,I)) ' time:' num2str(time(1))]);
% 
% save(['results\video_akiyo_fixtube_' num2str(kk) '.mat'])
% 
% end
% 
% 
% for kk=0.3:0.1:0.5
% 
% load(['results\video_foreman_randomstripe_' num2str(kk) '.mat'])
% disp(['I_OPTRC PSNR: ' num2str(psnr(I_OPTRC,I)) ' time:' num2str(time(1))]);
% option=get_option_default(MissM_int,Mask_int,'OPTRC',opt);
% t1=tic;
% I_OPTRC=OPTRC(MissM_int,Mask_int,I,option);
% I_OPTRC=(1-Mask).*I_OPTRC+Mask.*MissM;
% time(1)=toc(t1);
% disp(['I_OPTRC PSNR: ' num2str(psnr(I_OPTRC,I)) ' time:' num2str(time(1))]);
% 
% save(['results\video_foreman_randomstripe_' num2str(kk) '.mat'])
% 
% end
% 
% 
% for kk=0.3:0.1:0.5
% 
% load(['results\video_stefan_stripe_' num2str(kk) '.mat'])
% disp(['I_OPTRC PSNR: ' num2str(psnr(I_OPTRC,I)) ' time:' num2str(time(1))]);
% option=get_option_default(MissM_int,Mask_int,'OPTRC',opt);
% t1=tic;
% I_OPTRC=OPTRC(MissM_int,Mask_int,I,option);
% I_OPTRC=(1-Mask).*I_OPTRC+Mask.*MissM;
% time(1)=toc(t1);
% disp(['I_OPTRC PSNR: ' num2str(psnr(I_OPTRC,I)) ' time:' num2str(time(1))]);
% 
% save(['results\video_stefan_stripe_' num2str(kk) '.mat'])
% 
% end


