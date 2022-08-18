function X=TTNN(MissM,Mask,I,option)

beta=option.beta;
maxitr=option.maxitr;
stopc=option.stopc;
debug=option.debug;

%% initialize
X=MissM;
N=ndims(MissM);
Y=cell(N,1);
sk=[];
sz=size(MissM);
alpha=zeros(1,N-1);

for k=1:1:N-1
    M{k}=MissM;
    Y{k}=zeros(size(X));
    sk=[sk prod(sz(1:k))];
    alpha(k)=min(prod(sz(1:k)),prod(sz(k+1:N)));
end

alpha=alpha/sum(alpha);

for kk=1:1:maxitr
    
    X_pre=X;
    
    %% update M
    for k=1:1:N-1
        A=reshape(X,sk(k),[])+1/beta*reshape(Y{k},sk(k),[]); 
        [U,S,V]=svd(A,'econ');
        Sv=diag(S);
        Sv=max(Sv-1/beta,0);
        DA=U*diag(Sv)*V';
        M{k}=reshape(DA,sz);
    end
       
    %% update X
    X=0;
    for k=1:1:N-1
        X=X+alpha(k)*(1-Mask).*(M{k}-1\beta*Y{k});
    end
    X=X+Mask.*MissM;
    
    %% update Y
    for k=1:1:N-1
        Y{k}=Y{k}+beta*(X-M{k});
    end
    
    %beta=beta*1.05;

    E=X-X_pre;
    err=norm(E(:))/norm(MissM(:));
    
    if err<stopc
        break;
    end
    
    if debug&&mod(kk,10)==0
        err1=X-I;
        fprintf('Iter %.0f, Diff %.2f\n',kk,norm(err1(:)));
    end


end