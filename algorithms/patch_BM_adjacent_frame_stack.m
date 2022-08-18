function [pos_arr, error_arr] = patch_BM_adjacent_frame_stack(extractPatchnew, extractPatchMasknew, refpatch, refMask, totalpatchsize, param, patchind)

search_window_size = param.search_window_size; 
K=param.K;

[row,col]=ind2sub(totalpatchsize,patchind);
%// search window range <--> all searchable patch indices
ind_search_m = max(row-search_window_size, 1) : min(row+search_window_size, totalpatchsize(1));
ind_search_n = max(col-search_window_size, 1) : min(col+search_window_size, totalpatchsize(2));
[xxn,yyn] = meshgrid(ind_search_m,ind_search_n);
idx = sub2ind(totalpatchsize, xxn(:), yyn(:));

ind_ob=refMask==1;
refpatch=refpatch(ind_ob);
refMask=refMask(ind_ob);

newpatch = extractPatchnew(ind_ob,idx,:);
newMask = extractPatchMasknew(ind_ob,idx,:);

dis=single(refpatch)-single(newpatch);
coMask=single(refMask.*newMask);
metric=sum((coMask.*dis).^2,1)./sum(coMask);
metric=squeeze(metric);

[error, ind]=sort(metric); 

if size(ind,1)==1
    pos_arr  = idx(ind(1:K));
    error_arr  = error(1:K);
else
    pos_arr  = idx(ind(1:K,:));
    error_arr  = error(1:K,:);
end


end


