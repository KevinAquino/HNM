#!/bin/bash
#SBATCH --qos=vip
#SBATCH --job-name=OptimGenes
#SBATCH --mail-type=END
#SBATCH --mail-user=
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=2
#SBATCH --array=1-999
#SBATCH --output=Genes_A%A_%a.out
#SBATCH --error=Genes_A%A_%a.err

#Load Matlab 2017a module
ml MATLAB

matlab -nojvm -nodisplay<<-EOF

s=str2num(getenv('SLURM_ARRAY_TASK_ID'))
load APARC_genedata.mat
coefei=data;

ratio=coefei/(max(coefei)-min(coefei));
ratio=ratio-max(ratio)+1;

ratio(35:68)=ratio(1:34);

%% Read data..SC FC and time series of BOLD
%%%%%%%%%%  
load('SC_GenCog_PROB_30.mat');
load('DKatlas_timeseries.mat');

C=GrCV([1:34 42:75],[1:34 42:75]);
C=C/max(max(C))*0.2;

N=68;
NSUB=389;
NSUBSIM=389; 
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

FCtdata2=zeros(NSUB,N,N);
stau=zeros(NSUB,N);

%%%%%%%%%%%%%%
kk=1;
for nsub=1:NSUB
    BOLDdata=(squeeze(tsdata(:,:,nsub)))';
    for seed=1:N
        BOLDdata(seed,:)=BOLDdata(seed,:)-mean(BOLDdata(seed,:));
        timeseriedata(seed,:) =filtfilt(bfilt2,afilt2,BOLDdata(seed,:));
    end
    %% time decay
    for i=1:N
        for j=1:N
            FCtdata2(nsub,i,j)=corr2(timeseriedata(i,1:end-1)',timeseriedata(j,2:end)');
        end
    end
    for i=1:N
        ac=autocorr(timeseriedata(i,:));
        coef_lin_reg = polyfit(1:5,ac(4:8),1);
        tau(i) = -1/coef_lin_reg(1);
    end
    stau(nsub,:)=tau*delt;
    %%

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
FCtdata=squeeze(nanmean(FCtdata2,1));
staudata=nanmean(stau,1);


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

ALPHA=-0.9:0.05:0;
BETA=0:0.05:2;

[Ialpha Ibeta]=ind2sub([length(ALPHA),length(BETA)],s);
alpha=ALPHA(Ialpha);
beta=BETA(Ibeta);

PERTURB=0.:0.001:0.2;

SEED=[4 10 12 20];
sigfunc = @(A, x)(A(1) ./ (1 + exp(-A(2)*(x-A(3)))) + A(4));
options=optimset('MaxFunEvals',10000,'MaxIter',1000);

we=2.1;
gain=1+alpha+beta*ratio;
Jbal=Balance_J_gain(we,C,gain);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model FC and FCD
FCt2=zeros(NSUBSIM,N,N);
stausi=zeros(NSUBSIM,N);

kk=1;
for nsub=1:NSUBSIM
    neuro_act=zeros(round(1000*(Tmax+60)*0.754+1),N);
    sn=0.001*ones(N,1);
    sg=0.001*ones(N,1);
    nn=1;
    for t=0:dt:(1000*(Tmax+60)*0.754)
        xn=I0*Jexte+JN*w*sn+we*JN*C*sn-Jbal.*sg;
        xg=I0*Jexti+JN*sn-sg;
        rn=phie_gain(xn,gain);
        rg=phii_gain(xg,gain);
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
    
    bds=BOLD_act(20:754:end-10,:);
    FC_simul2(nsub,:,:)=corrcoef(bds);
    
    Tmax2=size(bds,1);
    Phase_BOLD_sim=zeros(N,Tmax2);
    BOLDsim=bds';
    for seed=1:N
        BOLDsim(seed,:)=BOLDsim(seed,:)-mean(BOLDsim(seed,:));
        signal_filt_sim =filtfilt(bfilt2,afilt2,BOLDsim(seed,:));
        timeserie(seed,:)=signal_filt_sim;
    end
        %% time decay
        
    for i=1:N
        for j=1:N
            FCt2(nsub,i,j)=corr2(timeserie(i,1:end-1)',timeserie(j,2:end)');
        end
    end
    for i=1:N
        ac=autocorr(timeserie(i,:));
        coef_lin_reg = polyfit(1:5,ac(4:8),1);
        tau(i) = -1/coef_lin_reg(1);
    end
    stausi(nsub,:)=tau*delt;
    %%

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
FCtsim=squeeze(nanmean(FCt2,1));
ErrDecay=corr2(diag(FCtdata),diag(FCtsim)); 
ErrDecay2=corr2(FCtdata(:),FCtsim(:)); 
stausim=nanmean(stausi,1);
ErrDecayBold=corr2(staudata,stausim);

FC_simul=squeeze(mean(FC_simul2,1));
cc=corrcoef(atanh(FC_emp(Isubdiag)),atanh(FC_simul(Isubdiag)));
FCfitt=cc(2); %% FC fitting

FCsim2=FC_simul-FC_simul.*eye(N);
GBCsim=mean(FCsim2,2);
GBCfitt=corr2(GBCemp,GBCsim);

[hh pp FCDfitt]=kstest2(cotsamplingdata,cotsamplingsim);  %% FCD fitting

%%%
save(sprintf('WgainPC_%03d.mat',s),'GBCfitt','FCfitt','FCDfitt');

EOF

