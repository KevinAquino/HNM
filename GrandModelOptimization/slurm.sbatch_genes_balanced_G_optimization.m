#!/bin/bash
#SBATCH --qos=vip
#SBATCH --job-name=OptimGenes
#SBATCH --mail-type=END
#SBATCH --mail-user=
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=2
#SBATCH --array=1-999
#SBATCH --output=Genes_A%A_%a.out
#SBATCH --error=Genes_A%A_%a.err

#Load Matlab 2017a module
ml MATLAB

matlab -nojvm -nodisplay<<-EOF

s=str2num(getenv('SLURM_ARRAY_TASK_ID'))

load DKcortex_selectedGenes.mat;
coefei=sum(expMeasures(:,18:25),2)./sum(expMeasures(:,2:6),2); %% ampa+nmda/gaba
ratioEI=coefei/(max(coefei)-min(coefei));
ratioEI=ratioEI-max(ratioEI)+1;

%ratioEI=(max(coefei)-coefei)./(max(coefei)-min(coefei));
%ratioEI=1-ratioEI;

ratioEI(35:68)=ratioEI(1:34);

%% Read data..SC FC and time series of BOLD
%%%%%%%%%%  
load('SC_GenCog_PROB_30.mat');
load('DKatlas_timeseries.mat');

C=GrCV([1:34 42:75],[1:34 42:75]);
C=C/max(max(C))*0.2;
N=68;
NSUB=389;
Tmax=616;
indexsub=1:NSUB;

for nsub=indexsub
    tsdata(:,:,nsub)=ts(:,[1:34 42:75],nsub);
    FCdata(nsub,:,:)=corrcoef(squeeze(tsdata(:,:,nsub)));
end

FC_emp=squeeze(mean(FCdata,1));

FCemp2=FC_emp-FC_emp.*eye(N);
GBCemp=mean(FCemp2,2);

Isubdiag = find(tril(ones(N),-1));

%%%%%%%%%%%%%%

flp = .008;           % lowpass frequency of filter
fhi = .08;           % highpass
delt = 0.754;            % sampling interval
k=2;                  % 2nd order butterworth filter
fnq=1/(2*delt);       % Nyquist frequency
Wn=[flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
[bfilt2,afilt2]=butter(k,Wn);   % construct the filter

%%%%%%%%%%%%%%
load('empiricalLEICAS_Melb.mat');
kk=1;
timeseriebolddata=[];
t_all=1;
for nsub=1:NSUB
    BOLDdata=(squeeze(tsdata(:,:,nsub)))';
    for seed=1:N
        BOLDdata(seed,:)=BOLDdata(seed,:)-mean(BOLDdata(seed,:));
        timeseriedata(seed,:) =filtfilt(bfilt2,afilt2,BOLDdata(seed,:));
 %       timeseriedataz(seed,:)=zscore(timeseriedata(seed,:));        
        Phase_BOLD_data(seed,:) = angle(hilbert(timeseriedata(seed,:)));
    end
 %   timeseriebolddata=horzcat(timeseriebolddata,timeseriedataz); 
    for t=1:Tmax
        iFC=zeros(N);
        for n=1:N
            for p=1:N
                iFC(n,p)=cos(Phase_BOLD_data(n,t)-Phase_BOLD_data(p,t));
            end
        end
        [V1,~]=eigs(iFC,1);
        Leading_Eig(t_all,:)=V1;
        Time_all(t_all)=s;
        t_all=t_all+1;
    end
    ii2=1;
    for t=1:18:Tmax-80
        jj2=1;
        cc=corrcoef((timeseriedata(:,t:t+80))');
        for t2=1:18:Tmax-80
            cc2=corrcoef((timeseriedata(:,t2:t2+80))');
            ca=corrcoef(cc(Isubdiag),cc2(Isubdiag));
            if jj2>ii2
                cotsamplingdata(kk)=ca(2);   %% this accumulate all elements of the FCD empirical
                kk=kk+1;
            end
            jj2=jj2+1;
        end
        ii2=ii2+1;
    end
end

Projectiondata=WHCP*Leading_Eig'; %timeseriebolddata;
Projectiondata(Projectiondata>3*std(Projectiondata,0,2))=1;
Projectiondata(Projectiondata<-3*std(Projectiondata,0,2))=-1;

for j=1:K
    distridata{j}=[];
end
for t=1:t_all-1
    for j=1:K
        di(j)=sqrt(sum((Leading_Eig(t,:)-CHCP(j,:)).^2)/N);
    end
    [aux indmin]=min(di);
    IDX(t)=indmin;
    distridata{indmin}=[distridata{indmin} aux];
end
for c=1:K
    Pdata(c)=mean(IDX==c);
end

%%%%%%%%%%%%%%%%%%
%% Parameters for the mean field model
dtt   = 1e-3;   % Sampling rate of simulated neuronal activity (seconds)
dt=0.1;

taon=100;
taog=10;
gamma=0.641;
sigma=0.01;
JN=0.15;
I0=0.382;
Jexte=1.;
Jexti=0.7;
w=1.4;

%%%%%%%%%%%%
%% Optimize
%%
WE=0:0.01:3;

we=WE(s);

PERTURB=0.:0.001:0.2;

SEED=[4 10 12 20];
sigfunc = @(A, x)(A(1) ./ (1 + exp(-A(2)*(x-A(3)))) + A(4));
options=optimset('MaxFunEvals',10000,'MaxIter',1000);

Jbal=Balance_J(we,C);

J=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model FC and FCD
kk=1;
timeseriebold=[];
t_all=1;
for nsub=1:NSUB
    neuro_act=zeros(round(1000*(Tmax+14)*0.754+1),N);
    sn=0.001*ones(N,1);
    sg=0.001*ones(N,1);
    nn=1;
    for t=0:dt:(1000*(Tmax+14)*0.754)
        xn=I0*Jexte+JN*w*sn+we*JN*C*sn-Jbal.*sg-J.*sg;
        xg=I0*Jexti+JN*sn-sg;
        rn=phie(xn);
        rg=phii(xg);
        sn=sn+dt*(-sn/taon+(1-sn)*gamma.*rn./1000.)+sqrt(dt)*sigma*randn(N,1);
        sn(sn>1) = 1;
        sn(sn<0) = 0;
        sg=sg+dt*(-sg/taog+rg./1000.)+sqrt(dt)*sigma*randn(N,1);
        sg(sg>1) = 1;
        sg(sg<0) = 0;
        j=j+1;
        if abs(mod(t,1))<0.01
            neuro_act(nn,:)=rn';
            nn=nn+1;
        end
    end
    nn=nn-1;
    
   %%%% BOLD empirical
    % Friston BALLOON MODEL
    T = nn*dtt; % Total time in seconds
    
    B = BOLD(T,neuro_act(1:nn,1)'); % B=BOLD activity, bf=Foutrier transform, f=frequency range)
    BOLD_act = zeros(length(B),N);
    BOLD_act(:,1) = B;
    
    for nnew=2:N
        B = BOLD(T,neuro_act(1:nn,nnew));
        BOLD_act(:,nnew) = B;
    end
    
    bds=BOLD_act(1:754:end,:);
    FC_simul2(nsub,:,:)=corrcoef(bds);
    
    Tmax2=size(bds,1);
    Phase_BOLD_sim=zeros(N,Tmax2);
    BOLDsim=bds';
    for seed=1:N
        BOLDsim(seed,:)=BOLDsim(seed,:)-mean(BOLDsim(seed,:));
        signal_filt_sim =filtfilt(bfilt2,afilt2,BOLDsim(seed,:));
        timeserie(seed,:)=signal_filt_sim;
    %    timeseriesim(seed,:)=zscore(signal_filt_sim);  
        Phase_BOLD_sim(seed,:) = angle(hilbert(signal_filt_sim));
    end
  %  timeseriebold=horzcat(timeseriebold,timeseriesim); 
    for t=1:Tmax2
       iFC=zeros(N);
       for n=1:N
           for p=1:N
               iFC(n,p)=cos(Phase_BOLD_sim(n,t)-Phase_BOLD_sim(p,t));
           end
       end
       [V1,~]=eigs(iFC,1);
       Leading_Eig(t_all,:)=V1;
       t_all=t_all+1;
    end   
    ii2=1;
    for t=1:18:Tmax2-80
        jj2=1;
        cc=corrcoef((timeserie(:,t:t+80))');
        for t2=1:18:Tmax2-80
            cc2=corrcoef((timeserie(:,t2:t2+80))');
            ca=corrcoef(cc(Isubdiag),cc2(Isubdiag));
            if jj2>ii2
                cotsamplingsim(kk)=ca(2);  %% FCD simulation
                kk=kk+1;
            end
            jj2=jj2+1;
        end
        ii2=ii2+1;
    end
end

FC_simul=squeeze(mean(FC_simul2,1));
cc=corrcoef(atanh(FC_emp(Isubdiag)),atanh(FC_simul(Isubdiag)));
FCfitt=cc(2); %% FC fitting

FCsim2=FC_simul-FC_simul.*eye(N);
GBCsim=mean(FCsim2,2);
GBCfitt=corr2(GBCemp,GBCsim);
    
Projectionsim=WHCP*(Leading_Eig(1:t_all-1,:))'; %timeseriebold;
Projectionsim(Projectionsim>3*std(Projectionsim,0,2))=1;
Projectionsim(Projectionsim<-3*std(Projectionsim,0,2))=-1;
for ass=1:NumAssemblies
   [aux paxu ksdist2(ass)]=kstest2(Projectiondata(ass,:),Projectionsim(ass,:));
end  
ICAfitt=mean(ksdist2);

[hh pp FCDfitt]=kstest2(cotsamplingdata,cotsamplingsim);  %% FCD fitting

for j=1:K
    distrisim{j}=[];
end
for t=1:t_all-1
    for j=1:K
        di(j)=sqrt(sum((Leading_Eig(t,:)-CHCP(j,:)).^2)/N);
    end
    [aux indmin]=min(di);
    IDX(t)=indmin;
    distrisim{indmin}=[distrisim{indmin} aux];
end
for c=1:K
    Psim(c)=mean(IDX==c);
end
    
for j=1:K
    if isempty(distrisim{j})
      distrisim{j}=zeros(1,1);
    end
end

for j=1:K
        [aux paxu ksdistclu2(j)]=kstest2(distridata{j},distrisim{j});
end
KSclu=sum(Pdata.*ksdistclu2.*(1+(sqrt(Pdata)-sqrt(Psim)).^2))/K;
    
KLPstates=sqrt(0.5*(sum((sqrt(Pdata)-sqrt(Psim)).^2)));

save(sprintf('WG_%03d.mat',s),'GBCfitt','FCfitt','FCDfitt');

EOF

