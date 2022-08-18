function Xhat=NLPTT_video_full(MissM,Mask,I,option)

Xhat=[];
fprintf('NLPTT Batch:     ');

if option.nChannel==3
    Xhat_temp=NLPTT(MissM(:,:,:,indp),Mask(:,:,:,indp),I,option);
else
    Xhat_temp=NLPTT(MissM(:,:,indp),Mask(:,:,indp),I,option);
end


end