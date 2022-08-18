function tensor_out = t_SVD_new(tensor, tensor_init, mask, para)

rho = 1;
if nargin < 4
    maxItr = 100;
else
    maxItr = para.maxItr;
end
eps = 1e-4;
mu = 10;

shape = size(tensor);

% x = zeros(n,1);
% z = zeros(n,1);
% u = zeros(n,1);

A = diag(sparse(double(~mask(:))));
b = A * tensor(:);
zInit = tensor_init(:);

z = zInit;
% uhat = zeros(n,1);
% init u
[U,~,V] = ntsvd(tensor_init,1);
E = zeros(size(U));
for i=1:size(U,3)
    E(:,:,i) = U(:,:,i)*V(:,:,i);
end
u = ifft(E,[],3); u = u(:);

P = double(~logical(diag(A)));
P = diag(sparse(P));
q = sparse(b);

for k = 1 : maxItr
    zold = z;

    x = P*(z - 1/rho*u) + q;     % x-update

    % z-update with relaxation
    [z, ~] = shrinkObj(x + 1/rho*u,1/rho,'tsvd_1',shape);
    u = u + rho*(x - z);
    
    r_norm = norm(x-z);
    s_norm = norm(-rho*(z - zold));
    
    if r_norm < eps && s_norm < eps
        break
    end

    Ek = (s_norm - r_norm)/mu;
    if Ek < -1
        rho = rho*2;
    elseif Ek > 1 && Ek < mu
        rho = rho/2;
    elseif Ek > mu
        disp('Pass!');
    end
end

tensor_out = reshape(x,shape);