function tensor_out = FBCP_inpaint(tensor, tensor_init, mask, para)
% mask --- is_missing_mat

if nargin < 4
    maxRank = 6;
else
    maxRank = para.maxRank;
end

maxItr = 500;

[model] = BayesCP_MP(tensor, 'obs', ~mask, 'init', 'rand', 'maxRank', maxRank, 'maxiters', maxItr, ...
	    'tol', 1e-5, 'dimRed', 1, 'verbose', 0, 'nd', 0.1);
tensor_out = double(ktensor(model.Z));   