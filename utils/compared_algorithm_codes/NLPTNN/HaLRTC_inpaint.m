function tensor_out = HaLRTC_inpaint(tensor, tensor_init, mask, para)
% mask --- is_missing_mat

if nargin < 4
    alpha = [1e-2, 1e-2, 1];
    maxItr = 200;
    rho = 1e-2;
else
    alpha = para.alpha;
    maxItr = 500;
    rho = 1e-5;
end


alpha = alpha / sum(alpha);

% parameters for t-SVD

[tensor_out, ~] = HaLRTC(...
    tensor, ...                       % a tensor whose elements in Observ are used for estimating missing value
    ~mask,...               % the index set indicating the obeserved elements
    alpha,...                  % the coefficient of the objective function, i.e., \|X\|_* := \alpha_i \|X_{i(i)}\|_*
    rho,...                      % the initial value of the parameter; it should be small enough
    maxItr,...               % the maximum iterations
    1e-6,...                 % the tolerance of the relative difference of outputs of two neighbor iterations
    tensor_init,...
    1);