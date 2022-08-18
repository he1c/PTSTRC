function M = T2M_r(T,i)
    % Input: Tensor T: k1 * 1 * k2
    % Output: Matrix M: k1 * k2
    T = TensPermute(T,i);
    M = reshape(T, size(T,1), []);
end