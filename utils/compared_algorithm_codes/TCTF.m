function TC = TCTF( MissM,Mask,I,option )

%% produce data
Omega=find(Mask(:)==1);
data=MissM(Omega);
known=Omega;
[n1,n2,n3]=size(MissM);
r_n=option.rank;

%% our method
if strcmp(option.method,'fixrank')
    option.rank_adj = 0*ones(1,n3);
    option.rank_min = r_n*ones(1,n3);
    EstCoreNway = r_n*ones(1,n3);
elseif strcmp(option.method,'image')
    option.rank_adj = 0*ones(1,n3);
    option.rank_min = [r_n,3*ones(1,n3-1)];
    EstCoreNway = [r_n,3*ones(1,n3-1)];
elseif strcmp(option.method,'video')
    option.rank_adj = 0*ones(1,n3);
    option.rank_min = [r_n,3*ones(1,5),3*ones(1,n3-1-5)];
    EstCoreNway = [r_n,3*ones(1,5),3*ones(1,n3-1-5)];%round(10*ones(1,n3));
else
    disp('error : unkown type of data')
    return;
end

Nway=[n1,n2,n3];
TC = TCTF_solver(data,known,I,Nway,EstCoreNway,option);

end

