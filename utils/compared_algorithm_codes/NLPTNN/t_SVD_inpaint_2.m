function tensor_out = t_SVD_inpaint_2(tensor, tensor_init, mask, para)
% mask --- is_missing_mat

alpha = 1.8;
if nargin < 4
    rho = 0.01;
    maxItr = 300;
else
    rho = para.rho;
    maxItr = para.maxItr;
end
myNorm = 'tSVD_1';
% parameters for t-SVD

dim = size(tensor);

A = diag(sparse(double(~mask(:))));
b = A * tensor(:);
bInit = zeros(size(b));
tensor_out = tensor_cpl_admm_init( A , b , bInit , rho , alpha , dim , maxItr , myNorm , 1); % QUIET = true
tensor_out = reshape(tensor_out,dim);