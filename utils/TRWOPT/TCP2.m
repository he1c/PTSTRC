function T2 = TCP2(U1,U2)
    % This is tensor connect product (TCP)
    % Input: U is a cell {U1, U2,..., Un}
    %        U_i is a tensor:       r_{i-1} * I_i * r_i
    % Output:T is a tensor:         r_0     * (I1...In) * r_n
    U{1}=U1;
    U{2}=U2;
    n = length(U);

    T = U{1};
    if(n==1)
        return;
    end
    
    for i =2 : n
        [lT, cT, ~] = size(T);
        [~, cU, rU] = size(U{i});
        T = reshape(left_unfold(T) * right_unfold(U{i}), [lT, cT*cU, rU]);
    end
    
    T2=reshape(permute(T,[2 3 1]),cT*cU,[]);
    
end