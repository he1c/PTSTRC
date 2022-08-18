% I=imread('data\124084.jpg');
% I=double(I)/255;
I=patch_full;

para_TR.max_tot     = 1e-6; para_TR.max_iter      =30;                     para_TR.disp =  1;
ObsRatio = 0.1;

P_Omega = zeros(size(I));
P_Omega( randsample( numel(I), round(ObsRatio * numel(I))) )=1; 
U_Omega = I.*P_Omega; 

MissM=patch_miss;
Mask=double(patch_miss~=-1);
MissM=MissM.*Mask;
nn=size(MissM);

N=ndims(MissM);

para_TR.r=4*ones(N,1);
tic
Utr_TR = Completion_TR_norm(MissM, Mask, para_TR);
toc
I_hat=fullTR(Utr_TR);
I_hat=reshape(I_hat,nn);

option.I=I;
option.r=4*ones(N,1);
tic
Utr_TR=TRSGD_scaled(MissM,Mask,option);
Ihat=fullTR(Utr_TR');
toc