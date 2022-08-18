function tensor_out = TV_inpaint(tensor, tensor_init, mask, para)

Observ = ~mask;
A=tensor;
[row, col, channel]=size(A);

index = cell(1,3);
for k = 1:channel
    [row_tem,col_tem] = find(Observ(:,:,k));
    index{k} = [row_tem';col_tem';linspace(k,k,length(row_tem))]';
end
index = [index{1}; index{2}; index{3}];
value = A(Observ);

tsize=[row, col, channel];
N=3;
alpha=[1/N, 1/N, 1/N];
beta=[1,1,0];

lambda_1=0.5;
lambda_2=1000;

tensor_out = LRTC_TV_II(index, value, lambda_1, lambda_2 ,alpha, beta, tsize, N );