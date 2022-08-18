function [Z1,AA,Ax] = Update_Z1_RLS(Mask, Z, MissM, AA, Ax, mu)
    
    n = length(Z);
    
    [rn, I1, r1] = size(Z{1});
    Z1 = zeros([rn, I1, r1]);

    B = TCP(Z(2:end)); 
    A = tens2mat(TensPermute(B,2), 1, [2,3]);
    P = tens2mat(Mask, 1, 2:n);  
    X = tens2mat(MissM, 1, 2:n);  
            
    for i=1:I1
        idx = P(i,:)==1;
        if sum(idx)>1
           AA{i}=mu*AA{i}+A(idx,:)'*A(idx,:);
           Ax{i}=mu*Ax{i}+A(idx,:)'*X(i,idx)';
        end
        Z1(:, i, :)=reshape(AA{i}\Ax{i}, [rn, 1, r1]);
    end
            
end 