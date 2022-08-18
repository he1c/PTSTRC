function tensor_out = t_SVD_121norm(tensor, tensor_init, mask, para)

rho1 = 1;
rho2 = 100;
if nargin < 4
    maxItr = 100;
else
    maxItr = para.maxItr;
end
eps = 1e-4;
mu = 10;

shape = size(tensor);

lambda=1/sqrt(shape(1));

X = tensor_init;

% uhat = zeros(n,1);
% init u
[T1,~,T2] = ntsvd(tensor_init,1);
E = zeros(size(T1));
for i=1:size(T1,3)
    E(:,:,i) = T1(:,:,i)*T2(:,:,i);
end
U = ifft(E,[],3);

[T1,~,T2] = ntsvd(tensor_init,1);
E = zeros(size(T1));
for i=1:size(T1,3)
    E(:,:,i) = T1(:,:,i)*T2(:,:,i);
end
V = ifft(E,[],3);

L=0;

S=0;

for k = 1 : maxItr
    Lold = L;
    Sold = S;
    
    % L-update with relaxation
    [L, ~] = shrinkObj(X(:) - 1/rho1*(-U(:) + S(:)),1/rho1,'tsvd_1',shape);
    L=reshape(L,shape);
    
    % update S
    Y = X+U-L;
    S=zeros(size(Y));
    for i=1:1:size(Y,1)
        for j=1:1:size(Y,3)
            temp=Y(i,:,j);
            S(i,:,j)=max(1-lambda/rho1/norm(temp(:)),0)*Y(i,:,j);
        end
    end
    
    % update Z
    Z = mask.*(X + 1/rho1*V) + ~mask.*tensor;     

    % update X
    X = (rho1*(L+S-U)+rho2*(Z-V))/(rho1+rho2);
    
    % update U,V
    U = U + rho1*(X-S-L);
    V = V + rho1*(X-Z);
      
    r_norm = norm(L(:)-Lold(:));
    s_norm = norm(S(:)-Sold(:));
    
    if r_norm < eps && s_norm < eps
        break
    end

    Ek = (s_norm - r_norm)/mu;
    if Ek < -1
        rho1 = rho1*2;
    elseif Ek > 1 && Ek < mu
        rho1 = rho1/2;
    elseif Ek > mu
        %disp('Pass!');
    end
end

tensor_out = X;