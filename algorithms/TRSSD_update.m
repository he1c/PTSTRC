function Z=TRSSD_update(MissM,Mask,Z_pre,newFramenum)

MissM=MissM.*Mask;
N=ndims(MissM);

% Initialization
Z = Z_pre;
t_order=ndims(Mask);
Z{end}=Z{end}(:,end-newFramenum+1:end,:);

%% avoid null patch
temp=MissM(:);
tempMask=Mask(:);
temp(tempMask==0)=[];
r_e=std(temp);
if r_e<1e-3
    return;
end

r=zeros(N,1);
for i=1:1:N
    r(i)=size(Z{i},1);
end
r(N+1)=r(1);

B = TCP(Z(1:end-1));              % r1 * (I2...In) * rn

for kk=1:1:newFramenum
    if t_order==3
        p = reshape(Mask(:,:,kk),[],1);  % I_1 * (I_2...I_n)
        x = reshape(MissM(:,:,kk),[],1); % I_1 * (I_2...I_n)   
    elseif t_order==4
        p = reshape(Mask(:,:,:,kk),[],1);  % I_1 * (I_2...I_n)
        x = reshape(MissM(:,:,:,kk),[],1); % I_1 * (I_2...I_n) 
    else
        disp('Wrong order!');
        return;
    end
    idx = find(p==1);
    if ~isempty(idx)
        A=tens2mat(TensPermute(B(:, idx, :),2), 1, [2,3]);
        b=x(idx);
        AA=A'*A;
        temp=(A'*A+1e-10*norm(AA(:))*eye(size(AA,2)))\(A'*b);
        Z{end}(:, kk, :) = reshape(temp, [r(end-1), 1, r(1)]);
    end
end


for itr=1:1:1

for k=1:1:N-1
   
    X=reshape(permute(Z{k},[2 1 3]),size(Z{k},2),[]);   
    
    nk=[k+1:N 1:k-1];
    Znk=TCP(Z(nk));
    Y=reshape(TensPermute(Znk,2),size(Znk,2),[])';

    MissMk=reshape(TensPermute(MissM, k),size(MissM,k),[]);
    Maskk=reshape(TensPermute(Mask, k),size(Mask,k),[]);
    
    dY=-(MissMk-Maskk.*(X*Y))*Y';
    YY=Y*Y';
    dY2=-dY/(YY+1e-10*norm(YY(:))*eye(size(YY,2)));
    dY3=Maskk.*(dY2*Y);  
    tx=trace(dY'*dY2)/norm(dY3(:)).^2;
    X=X-tx*dY2;

    Z{k}=permute(reshape(X,[],r(k),r(k+1)),[2 1 3]);
    
end


end

end