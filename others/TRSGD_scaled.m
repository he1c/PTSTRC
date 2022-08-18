function Z=TRSGD_scaled(MissM,Mask,option)

r=option.r;
MissM=MissM.*Mask;

N=ndims(MissM);
r(N+1)=r(1);

Z=cell(1,N);
for k=1:1:N
    Z{k}=0.1*rand(r(k),size(MissM,k),r(k+1));
end
%Z = TR_Initialization0(MissM, r)';

for itr=1:1:30

for k=1:1:N
   
    X=reshape(permute(Z{k},[2 1 3]),size(Z{k},2),[]);   
    
    nk=[k+1:N 1:k-1];
    Znk=TCP(Z(nk));
    Y=reshape(TensPermute(Znk,2),size(Znk,2),[])';

    MissMk=reshape(TensPermute(MissM, k),size(MissM,k),[]);
    Maskk=reshape(TensPermute(Mask, k),size(Mask,k),[]);
    
    dY=-(MissMk-Maskk.*(X*Y))*Y';
    dY2=-dY*inv(Y*Y');
    dY3=Maskk.*(dY2*Y);  
    tx=trace(dY'*dY2)/norm(dY3(:)).^2;

%      dY=-(MissMk-Maskk.*(X*Y))*Y';
%      dY2=Maskk.*(dY*Y);
%      tx=trace(dY'*dY)/norm(dY2(:)).^2;

    X=X-tx*dY2;
    Z{k}=permute(reshape(X,[],r(k),r(k+1)),[2 1 3]);

%  normalization
%     if k<N
%         U_temp=reshape(TensPermute(Z{k},2),size(Z{k},2),[]);
%         for j=1:1:size(U_temp,2)
%             U_temp(:,j)=U_temp(:,j)./norm(U_temp(:,j));
%         end
%             Z{k}=TensPermute(reshape(U_temp,size(Z{k},2),size(Z{k},1),size(Z{k},3)),3);
%     end


end

% Ihat=fullTR(Z');
% disp(psnr(Ihat,I));

end

end