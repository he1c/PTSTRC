function [Xhat,ind_all]=NLPTT_video_10(MissM,Mask,I,option)

if option.framenum<10
    Xhat=NLPTT(MissM,Mask,I,option);
else
    Xhat=[];
    ind_all=[];
    ind=1:1:option.framenum_all;
    fprintf('NLPTT Batch:     ');
    for i=1:3:option.framenum_all/10%length(ind)
        %fprintf('\b\b\b\b\b%5i',i);
        indp=(i-1)*10+1:min(i*10,option.framenum_all)
        if option.nChannel==3
            Xhat_temp=NLPTT(MissM(:,:,:,indp),Mask(:,:,:,indp),I,option);
        else
            Xhat_temp=NLPTT(MissM(:,:,indp),Mask(:,:,indp),I,option);
        end
        Xhat=cat(ndims(MissM),Xhat,Xhat_temp);
        ind_all=[ind_all indp];
    end
    fprintf('\n');
end

end