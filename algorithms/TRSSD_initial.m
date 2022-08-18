function Z=TRSSD_initial(MissM,Mask,option)

MissM=MissM.*Mask;
N=ndims(MissM);
Z=cell(N,1);

if option.adaptiverank==1
    %% adaptive rank estimation
    temp=MissM(:);
    tempMask=Mask(:);
    p=sum(tempMask(:)~=0)/length(tempMask(:));
    temp(tempMask==0)=[];
    r_e=var(temp);
    m_e=mean(temp)^(1/N);
    if r_e<1e-6        
        for k=1:1:N
            Z{k}=m_e*ones(1,size(MissM,k),1);
        end
        return;
    end
    r_a=min(ceil(r_e*sqrt(p)*1000),floor(sqrt(p)*6)+4);%min(round(p*10)+4,option.r(1)));
    r=r_a*ones(N,1);
else
    r=option.r;
end

r(N+1)=r(1);
for k=1:1:N
    Z{k}=m_e/max(r)+0.1*randn(r(k),size(MissM,k),r(k+1));
end

M_pre=0;

for itr=1:1:option.max_iter

for k=1:1:N
   
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

    M=X*Y;
    
    err=norm(M(:)-M_pre(:))/norm(M_pre(:));
    
    if err<option.stopc
        break;
    end
    
    M_pre=M;
    

end

end