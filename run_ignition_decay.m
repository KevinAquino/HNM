% This code runs the ignition and the decay variables
function run_ignition_decay(slurm_id,model_variant,range_perturb,jbal_flag)

% Default is to run the code to find the optimal Jbalanced terms
if(nargin<4)
	jbal_flag=0;
end

STEPS_P=4; % this is how many steps you skip
% Default is to have a full range of the perturbation, otherwise break it up!
if(nargin<3 || isempty(range_perturb))
	PERTURB=0.:0.001:0.2;
else
	PERTURB=0.:0.001:0.2;

	% Steps of 4
	range_perturb=str2num(range_perturb);
	PERTURB=PERTURB(range_perturb:STEPS_P:end);


	% And a little extra for the first one to include 0.2 stimulation
	if(range_perturb==1)
		PERTURB=[PERTURB,0.2];
	end
end

run_for_fig=1;
if(run_for_fig)
	PERTURB=0.2;
end
% Now set up the Random number generator so we can reproduce the noisy terms:
% Also convert from a str to a num
slurm_id=str2num(slurm_id);
rng(slurm_id);

% Defining the SEED
SEED=10;

% First load in all the variables again

switch model_variant
	case 'EI'		
		% Load in the genes 
		load DKcortex_selectedGenes.mat;		
		coefe=sum(expMeasures(:,18:25),2);
		ratioE=coefe/(max(coefe));
		ratioE(35:68)=ratioE(1:34);

		coefi=sum(expMeasures(:,[2:9 12:14]),2); %% 18:21 ampa+ 22:25 nmda/gaba
		ratioI=coefi/(max(coefi));
		ratioI(35:68)=ratioI(1:34);
		ratio=ratioE./ratioI;
		ratio=ratio/(max(ratio)-min(ratio));
		ratio=ratio-max(ratio)+1;

		% Grand model simulation parameters
		alpha=-0.3;  %% -0.3
		beta=1.8;   %%1.8
		we=2.1;

	case 'T1T2'
		load myelin_HCP_dk68.mat;
		ratio=t1t2Cortex';
		ratio=ratio/(max(ratio)-min(ratio));
		ratio=ratio-max(ratio)+1;
		ratio(find(ratio<0))=0;
		alpha=-0.7;  %% -0.3
		beta=1.4;   %%1.8
		we=2.1;

	case 'PCA'
		% Here is the PCA model
		load APARC_genedata.mat
		coefei=data;

		ratio=coefei/(max(coefei)-min(coefei));
		ratio=ratio-max(ratio)+1;

		ratio(35:68)=ratio(1:34);
		alpha=-0.75;  %% -0.3
		beta=1.6;   %%1.8
		we=2.1;

	case 'Null'
		% Loads in the null spatial model for the genes		
		% indrnd=randi(10000);
		% Loading in the pseudo random vector, to have more control over this and to be able to pre-calculate the values
		load null_pseudoRandom.mat;

		indrnd=randVector(slurm_id);

		load DKcortex_selectedGenes_Permuted.mat;		
		% Load one of the permuted genes from a spin test:
		coefe=sum(expMeasuresPerm(:,18:25,indrnd),2);
		ratioE=coefe/(max(coefe));

		% Reflect the left to the right hemisphere
		ratioE(35:68)=ratioE(1:34);
		% Here perform normalization
		coefi=sum(expMeasuresPerm(:,[2:9 12:14],indrnd),2); %% 18:21 ampa+ 22:25 nmda/gaba
		ratioI=coefi/(max(coefi));
		% Also refelect this measure
		ratioI(35:68)=ratioI(1:34);
		% Now calcualte the ratio between excitatory genes / inhibitory genes
		ratio=ratioE./ratioI;

		% Perform some ratio normalization:
		ratio=ratio/(max(ratio)-min(ratio));
		ratio=ratio-max(ratio)+1;

		% Last of all the parameters for the simulation - derived from a different simulation that performs and exhaustive search over these params
		we=2.1;
		alpha=-0.35;
		beta=1.85;		

	case 'BEI'
		ratio=ones(68,1);
		alpha=0;
		beta=0;	
		we=2.1;
	end

	
	% Define now the gain for the models:
	gain=1+alpha+beta*ratio;

	% Define the SC matrix
	%%%%%%%%%%  
	load('SC_GenCog_PROB_30.mat');

	C=GrCV([1:34 42:75],[1:34 42:75]);
	C=C/max(max(C))*0.2;


	% Balancing the model - re-calculating or re-loading, i have a parameter that you can switch it to run it
	% or preload it - for the moment, if you want to calculate it, it saves jbal then exists
	if(jbal_flag==1)
		% Run the Jbalance optimization
		Jbal=Balance_J_gain(we,C,gain);
		system(['mkdir -p ',model_variant])
		switch model_variant
			case 'Null',
				save([model_variant,'/JBAL_',model_variant,'_',num2str(slurm_id),'.mat'],'Jbal');				
			otherwise
				save([model_variant,'/JBAL_',model_variant,'.mat'],'Jbal');
		end
		return
	else
		switch model_variant
			case 'Null',
				load([model_variant,'/JBAL_',model_variant,'_',num2str(slurm_id),'.mat']);
			otherwise
				load([model_variant,'/JBAL_',model_variant,'.mat']);
		end

	end

		
	% The filename to save the simulation
	if(strcmp(computer,'MACI64'))
		cluster_folder='./';
	else
		% Please put in the path you need
		disp('Enter your cluster folder on this line... exiting..')		
		return
		% cluster_folder='/myClusterFolder/';
	end
    
	filename=[cluster_folder,model_variant,'/',model_variant,'_',num2str(slurm_id),'_Peturb_range_',num2str(range_perturb)];

	% Here a little flag - if it was already run, do not run it again (can change later)
	if(exist([filename,'.mat']))
		disp('!!!!!!!!Found simulation - NOT overriding -- exiting now!')
		return
	end

	%%%%%%%%%%%%%%%%%%
	% Parameters for the mean field model:
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



	N=68; 
	Tmax=616;


	% Settings for the ignition code:
	sigfunc = @(A, x)(A(1) ./ (1 + exp(-A(2)*(x-A(3)))) + A(4));
	options=optimset('MaxFunEvals',10000,'MaxIter',1000);

	nseed=1;
	Tmax=7000;
	SUBTR=500;
	nump=length(PERTURB);
	neuro_act2=zeros(SUBTR,Tmax+1,N);



	for seed=SEED   %% Only 1 seed is being stimulated but here for convenience.	    	   
	    tseed=1:N;
	    tseed(find(seed==tseed))=[];
	    tseed(find(N/2+seed==tseed))=[];
	    kk=1;
	    for perturb=PERTURB  %% for each stimulation..different strengths (Istim)	    	
	        Istim=zeros(N,1);


	        % Here peforming the simulation over SUBTR times to form an average - this forms all the information for the next sections

	        for sub=1:SUBTR
	        	disp(['Running ... ',num2str(sub)]);
	        	% Initalize the code and perform the stimulation:
	            nn=1;
	            sn=0.001*ones(N,1);
	            sg=0.001*ones(N,1);
	            for t=0:dt:Tmax
	                if t==3000
	                    Istim(seed)=perturb;
	                    Istim(N/2+seed)=perturb;
	                elseif t==4000
	                    Istim(seed)=0;
	                    Istim(N/2+seed)=0;
	                end
	                xn=Istim+I0*Jexte+JN*w*sn+we*JN*C*sn-Jbal.*sg;
	                xg=I0*Jexti+JN*sn-sg;
	                rn=phie_gain(xn,gain);
	                rg=phii_gain(xg,gain);
	                sn=sn+dt*(-sn/taon+(1-sn)*gamma.*rn./1000.)+sqrt(dt)*sigma*randn(N,1);
	                sn(sn>1) = 1;
	                sn(sn<0) = 0;
	                sg=sg+dt*(-sg/taog+rg./1000.)+sqrt(dt)*sigma*randn(N,1);
	                sg(sg>1) = 1;
	                sg(sg<0) = 0;
	                if abs(mod(t,1))<0.01
	                    neuro_act2(sub,nn,:)=rn';
	                    nn=nn+1;
	                end
	            end
	        end


	        % From all the sub THR, we now look at averaging the activity so that it is all saved!

	        neuro_act1=squeeze(mean(neuro_act2(1:SUBTR,:,:),1));
	        ntwi=1;
	        for twi=1:20:nn-1-19
	            neuro_actf(ntwi,:)=mean(neuro_act1(twi:twi+19,:));
	            ntwi=ntwi+1;
	        end



	        ssnum=1;
	        for ss=tseed
	            peakrate3(ssnum)=mean(neuro_actf(160:200,ss),1);
	            basal(ssnum)=mean(neuro_actf(100:140,ss));
	            decayneuro=squeeze(neuro_actf(200:end,ss));
	        end
	        peakrate2(kk,:)=peakrate3./basal;
	        neuro_act_all(:,:,kk) = neuro_actf;
	        kk=kk+1;
	    end
	    
	    % Now save the results, but make a folder (doesnt overwrite it)
	    system(['mkdir -p ',cluster_folder,model_variant])
	    
	    save(filename,'neuro_act_all','PERTURB','model_variant','slurm_id');


	end





end

