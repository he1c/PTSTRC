function Z=TRSGD_scaled_add(MissM,Mask,option,Z_pre,newFramenum, maxsize)

r=option.r;
MissM=MissM.*Mask;

N=ndims(MissM);
r(N+1)=r(1);

% Initialization
Z = Z_pre;
n_f_old=size(Z{end},2);
t_order=ndims(Mask);

B = TCP(Z(1:end-1));              % r1 * (I2...In) * rn

for kk=1:1:newFramenum
    if t_order==3
        p = reshape(Mask(:,:,end-newFramenum+kk),[],1);  % I_1 * (I_2...I_n)
        x = reshape(MissM(:,:,end-newFramenum+kk),[],1); % I_1 * (I_2...I_n)   
    elseif t_order==4
        p = reshape(Mask(:,:,:,end-newFramenum+kk),[],1);  % I_1 * (I_2...I_n)
        x = reshape(MissM(:,:,:,end-newFramenum+kk),[],1); % I_1 * (I_2...I_n) 
    else
        disp('Wrong order!');
        return;
    end
    idx = find(p==1);
    Z{end}(:, n_f_old+kk, :) = reshape(tens2mat(TensPermute(B(:, idx, :),2), 1, [2,3])\x(idx), [option.r(end), 1, option.r(1)]);
end
if size(Z{end},2)>maxsize
    Z{end}=Z{end}(:,newFramenum+1:end,:);
    Z_pre{end}=Z_pre{end}(:,newFramenum+1:end,:);
end

MissM=MissM(:,:,:,end-newFramenum+1:end);

M_pre=reshape(T2M_r(B,2)*T2M_n(Z_pre{end})',[size(MissM,1) size(MissM,2) size(MissM,3) size(Z_pre{end},2)]);

MissM=cat(4,M_pre,MissM);

for k=1:1:N-1
   
    X=reshape(permute(Z{k},[2 1 3]),size(Z{k},2),[]);   
    
    nk=[k+1:N 1:k-1];
    Znk=TCP(Z(nk));
    Y=reshape(TensPermute(Znk,2),size(Znk,2),[]);

    MissMk=reshape(TensPermute(MissM, k),size(MissM,k),[]);
    Maskk=reshape(TensPermute(Mask, k),size(Mask,k),[]);
    
    dY=-(Maskk.*MissMk-Maskk.*(X*Y'))*Y;
    dY2=-dY*inv(Y'*Y);
    dY3=Maskk.*(dY2*Y');  
    tx=trace(dY'*dY2)/norm(dY3(:)).^2;

    X=X-tx*dY2;
    Z{k}=permute(reshape(X,[],r(k),r(k+1)),[2 1 3]);

end

end