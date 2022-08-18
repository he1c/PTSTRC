function tensor_out = NLR_inpaint(tensor, tensor_init, mask, para)

% [h,w] = size(tensor);
% rate = sum(~mask(:))/numel(mask);
% par = Set_parameters(rate, rate*h, 1);
% par.s_model = 1;
% 
% P = find(mask);
% A = @(z) A_fhp(z, P, h, w);
% At = @(z) At_fhp(z, P, h, w);
% 
% par.y = A(tensor(:));
% par.ori_im = tensor;
% par.picks = P;
% 
% [tensor_out,~,~] = NLR_CS_Reconstruction( par, A, At ); 

tensor_out = Image_CS( tensor, 1, sum(~mask(:))/numel(mask), 0 );

tensor_out = real(tensor_out);
tensor_out(~mask) = tensor(~mask);