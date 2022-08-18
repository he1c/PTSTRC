function U1 = Updata_U1_CP(Omega, U, X_Omega)
    % Input
    % Omega:    Observe idx  I1 * ... * In
    % X_Omega:  Missing data I1 * ... * In
    % U:        {U1,..., Un}, each Ui is a 3 mode tensor
    %
    % Output
    % U1:       update U1 to solve:
    % U1 =      argmin ||Omega \circ f(yU2..Ud) - X_Omega||_F^2

    A = kr(U(2:end));              % r1 * (I2...In) * rn
    
    P = reshape(Omega,size(Omega,1),[]);  % I_1 * (I_2...I_n)
    X = reshape(X_Omega,size(X_Omega,1),[]);  % I_1 * (I_2...I_n)
    
    [I1, r1] = size(U{1});
    U1 = zeros([I1, r1]);
    
    for i=1:I1
        idx = P(i,:)==1;
        if sum(idx)>1
            U1(i, :) = reshape(A(idx,:)\(X(i, idx))', [1, r1]);
        end
    end
        
end 