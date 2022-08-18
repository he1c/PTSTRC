function Z1 = Update_Z1_full(Mask, Z, MissM)

    B = TCP(Z(2:end));             
    
    n = length(Z);
    P = tens2mat(Mask,   1, 2:n);  
    X = tens2mat(MissM, 1, 2:n);  
    
    [rn, I1, r1] = size(Z{1});
    Z1 = zeros([rn, I1, r1]);
         
    A=tens2mat(TensPermute(B,2), 1, [2,3]);
    
    parfor i=1:I1
        idx = P(i,:)==1;
        if sum(idx)>1
            Z1(:, i, :) = reshape(A(idx,:)\(X(i, idx))', [rn, 1, r1]);
        end
    end
        
end 