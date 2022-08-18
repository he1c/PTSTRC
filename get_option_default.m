function option=get_option_default(MissM,Mask,alg,opt)

option=[];

framenum=opt.framenum;

if ndims(MissM)==3&&framenum~=1
    option.nChannel=1;
else
    option.nChannel=3;
end

option.stopc=1e-3;
option.maxitr  = 200;
option.debug=0;

switch alg
        
     case 'OPTRC'      
        option.border_ext=20;
        option.interval=3;
        option.dim=36;
        option.search_window_size=20;     
        option.stride=24; 
        option.se=strel('square',3);
        option.buffersize=2;
        option.p_frame=30;
        option.p_frame_o=10;
        option.TRrank = 8;
        option.TRmaxiter = 10;
        option.TRstopc  = 1e-2;
        option.TRdisp =  0;
        option.adaptiverank = 1;
        option.min_ob_num=5;
        option.overlap=3;
        option.mistrack=1600;
        
    case 'PTRC'
        n_t=opt.n_t;
        N_t=length(n_t);
        L=ceil(N_t/2);
        p=sum(double(Mask(:)==1))/sum(double(Mask(:)));
        sk=[];
        for nn=1:N_t
            order=[nn:N_t 1:nn-1];
            M=reshape(MissM,prod(n_t(order(1:L))),[]);
            sk=[sk max(ceil((min(size(M)))*opt.eta*sqrt(p)),max(floor(sqrt(framenum)*2),3))];
        end
        option.d=ceil(N_t/2);
        option.beta=1/N_t*ones(N_t,1);
        option.r=sk;
        option.lambda=1;
        
    case 'TRNN'
        option.beta = 1e-5;
        option.alpha = 1e-4;
        
    case 'TNN'
        option.rho=0.01;
        option.alpha=1;
        option.size=size(MissM);
        option.myNorm= 'tSVD_1';
        
    case 'TCTF'
        if framenum==1
            option.rank   = 10;
            option.method   = 'image';
        else
            option.rank   = 20;
            option.method   = 'video';
        end
        option.stopc=1e-5; 
        
    case 'TRLRF'
        option.r=5*ones(1,ndims(MissM)); % TR-rank 
        option.maxitr=300; % maxiter 300~500
        option.tol=1e-6; % 1e-6~1e-8
        option.Lambda=5; % usually 1~10
        option.ro=1.1; % 1~1.5
        option.K=1e0; % 1e-1~1e0 
        
    case 'TMAC'
        option.maxit = 300; 
        option.tol = 1e-5; % run to maxit by using negative tolerance
        option.Mtr = MissM; % pass the true tensor to calculate the fitting
        option.alpha_adj = 0;
        option.rank_adj = -1*ones(1,3);
        if framenum==1
            option.rank_min = 10*ones(1,3);
            option.rank_max = 10*ones(1,3);
        else
            option.rank_min = [100,100,50];%60*ones(1,3);
            option.rank_max = [100,100,50];%60*ones(1,3);
        end
            
    case 'OLSTEC'
        option.maxepochs       = 2;
        option.tolcost         =  1e-8;
        option.lambda          = 0.01;
        option.stepsize        = 0.1;
        option.mu              = 0.05;
        option.permute_on      = false;    
        option.store_subinfo   = true;     
        option.store_matrix    = true; 
        option.verbose         = 0; 
        option.rank  = 100;
        option.size = size(MissM);
        option.Gamma_in=[];
        option.X_init=[];
        
    case 'NLPTNN'
        option.radius = 18;
        option.search_rad = 20;
        option.search_gap = 1;
        option.THRESHOLD = 19;
        
    case 'NLPTT'
        option.K=20;
        option.dim=36;
        option.stride=24;  
        option.beta=0.03;
        option.maxitr=100;
        option.stopc=1e-3;
        option.framenum=framenum;

end


end