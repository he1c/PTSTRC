function Utr = Completion_TR_test(U_Omega, P_Omega, para,I)
    % U_Omega: I1 * I2 *... *In tensor with missing entries filled by 0
    % P_Omega: I1 * I2 *... *In binary observation tensor
    % para:    para.max_iter, para.max_tot, para.r
    %
    % Utr:     Decomposition term
    
    % Initialization
    Utr = TR_Initialization0(U_Omega,  para.r);
    dummy   = Utr{end};
    
    d = length(para.r);    
    tot     =1;
    iter    =0;
    
    while(tot>=para.max_tot && iter <=para.max_iter)
        for i=1:d
            Utr{i} = Updata_U1_temp(TensPermute(P_Omega, i), Utr([i:d, 1:i-1]),TensPermute(U_Omega, i));
        end
        % Change of the last term as error
        tot = norm(T2V(Utr{d} - dummy))/norm(T2V(dummy));
        dummy = Utr{end};
        iter =iter+1;
        if para.disp ==1
            disp(['At iteration ', num2str(iter), ', the last term change is ', num2str(tot)]);
            TCP2=TCP(Utr(2:end));
            Ihat=reshape(T2M_n(Utr{1})*T2M_r(TCP2,2)',size(I));
            disp(psnr(Ihat,I));
        end
    end
end