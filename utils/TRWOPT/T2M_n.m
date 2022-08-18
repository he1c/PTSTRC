function M = T2M_n(T)
    % Input: Tensor T: k1 * 1 * k2
    % Output: Matrix M: k1 * k2
    T = permute(T,[2 1 3]);
    M = reshape(T,size(T,1),[]);
end