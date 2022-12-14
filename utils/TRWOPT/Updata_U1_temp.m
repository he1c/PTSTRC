function U1 = Updata_U1_temp(Omega, U, X_Omega)
    % Input
    % Omega:    Observe idx  I1 * ... * In
    % X_Omega:  Missing data I1 * ... * In
    % U:        {U1,..., Un}, each Ui is a 3 mode tensor
    %
    % Output
    % U1:       update U1 to solve:
    % U1 =      argmin ||Omega \circ f(yU2..Ud) - X_Omega||_F^2

    B = TCP(U(2:end));              % r1 * (I2...In) * rn
    %B = Z_neq(U,1);
    
    n = length(U);
    P = tens2mat(Omega,   1, 2:n);  % I_1 * (I_2...I_n)
    X = tens2mat(X_Omega, 1, 2:n);  % I_1 * (I_2...I_n)
    
    [rn, I1, r1] = size(U{1});
    if ismatrix(U{1})
       [rn, I1, r1] = size(reshape(U{1},1,[]));
    end
    U1 = zeros([rn, I1, r1]);
         
    A=tens2mat(TensPermute(B,2), 1, [2,3]);
    
    for i=1:I1
        idx = P(i,:)==1;
        if sum(idx)>1
            U1(:, i, :) = reshape(A(idx,:)\(X(i, idx))', [rn, 1, r1]);
        end
    end
        
end 