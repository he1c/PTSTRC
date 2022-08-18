

I=imread('D:\Postdoc\MATRC\flamingo\00000.jpg');
I=double(I)/255;

para_TR.max_tot     = 1e-3; para_TR.max_iter      =100;                     para_TR.disp =  1;
ObsRatio = 0.1;

P_Omega = zeros(size(I));
P_Omega( randsample( numel(I), round(ObsRatio * numel(I))) )=1; 
U_Omega = I.*P_Omega; 

idx=1:900;
MissM=U_Omega(idx,idx,:);
Mask=P_Omega(idx,idx,:);
nn=size(MissM);

n_t=[30,30,30,30,3];
MissM_t=reshape(MissM,n_t);
Mask_t=reshape(Mask,n_t);
N=ndims(MissM_t);

para_TR.r=20*ones(N,1);
tic
Utr_TR = Completion_TR(MissM_t, Mask_t, para_TR);
time_TR=toc;
I_hat=fullTR(Utr_TR);
I_hat=reshape(I_hat,nn);


%% ==============PTRC+TSC================
sk=[];
d=ceil(N/2);
for k=1:N
    order=[k:N 1:k-1];
    M_temp=reshape(MissM_t,prod(n_t(order(1:d))),[]);
    sk=[sk max(ceil((min(size(M_temp)))*0.4*sqrt(ObsRatio)),max(floor(sqrt(n_t(end))*2),3))];
end
optiontc.d = d;
optiontc.beta = 1/N*ones(N,1);
optiontc.r = sk;
optiontc.lambda = 1;
optiontc.stopc = 1e-3;
optiontc.maxitr  = 300;
optiontc.debug= 0 ;
I_hat2=PTRC(MissM_t,Mask_t,[],optiontc);
I_hat2=reshape(I_hat2,nn);


%% ==============PTRC+TSC================;
optiontc.r = 50*ones(N,1);
I_hat3=PTRC(MissM_t,Mask_t,[],optiontc);
I_hat3=reshape(I_hat3,nn);
% 
% % generate init data
% rank=10;
% Xinit.A = randn(n_t(1), rank);
% Xinit.B = randn(n_t(2), rank);    
% Xinit.C = randn(n_t(3), rank); 
% 
% 
% %% CPOPT (batch)
% clear options;   
% options.maxepochs       = 1;
% options.display_iters   = 1;
% options.store_subinfo   = true;     
% options.store_matrix    = true; 
% options.verbose         = 2; 
% 
% [Xsol_cp_wopt, ~, ~] = cp_wopt_mod(MissM_t, Mask_t, [], n_t, rank, []);
% I_hat4=Xsol_cp_wopt.data;















