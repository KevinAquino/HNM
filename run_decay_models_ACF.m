function run_decay_models_ACF(slurm_id,model_variant,jbal_flag,SEED)
    
    % Default is to run the code to find the optimal Jbalanced terms
    if(nargin<4)
        jbal_flag=0;
    end
    %
    % Now have the right slurm ID to start the random number generator
    slurm_id=str2num(slurm_id);
    rng(slurm_id);
    SEED=str2num(SEED);


    % Now for the five differet models set up the initial conditions
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
            sigma=1e-8;

        case 'T1T2'
            load myelin_HCP_dk68.mat;
            ratio=t1t2Cortex';
            ratio=ratio/(max(ratio)-min(ratio));
            ratio=ratio-max(ratio)+1;
            ratio(find(ratio<0))=0;
            alpha=-0.7;  %% -0.3
            beta=1.4;   %%1.8
            we=2.1;
            sigma=1e-8;

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
            sigma=1e-5;

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
            sigma=1e-8;

        case 'BEI'
            ratio=ones(68,1);
            alpha=0;
            beta=0; 
            we=2.1;
            sigma=1e-4;
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
        cluster_folder='/scratch/kg98/kevo/Ignition_decay/';
        filename=[cluster_folder,model_variant,'/',model_variant,'_',num2str(slurm_id),'_ACF_DECAY_',num2str(SEED),];

        % Here a little flag - if it was already run, do not run it again (can change later)
%         if(exist([filename,'.mat']))
%             disp('!!!!!!!!Found simulation - NOT overriding -- exiting now!')
%             return
%         end

        %%%%%%%%%%%%%%%%%%
        % Parameters for the mean field model:
        dtt   = 1e-3;   % Sampling rate of simulated neuronal activity (seconds)
        dt=0.1;
        taon=100;
        taog=10;
        gamma=0.641;
%         sigma=0.001;
%         sigma=1e-6;
        JN=0.15;
        I0=0.382;
        Jexte=1.;
        Jexti=0.7;
        w=1.4;



        N=68; 


        %% Check the temporal decay due to stimulation of the occipital cortex.
        % Parameters:
        nseed=1; % Seed counter
        Tmax=50000; % Maximal stimulation time (its a stimulation)
        SUBTR=30; % Stimulation iterations - so its run 30 times
        PERTURB=0.2; % Stimulation
        neuro_act2=zeros(SUBTR,Tmax+1,N);





    % ACF DECAY: First take the SEED region and place noise on them and simulate (same as Wang et al.)

        for seed=SEED   %% check all single area stimulation
            seed
            clear peakrate1 peak1 peakrate3 peakrate2 basal tdecay ignition1;
            tseed=1:N;
            kk=1;
            for perturb=PERTURB(end)  %% for each stimulation..different strengths (Istim)
                perturb
                Istim=zeros(N,1);
                for sub=1:SUBTR
                    disp(['Simulation: ',num2str(sub)]);
                    nn=1;
                    sn=0.001*ones(N,1);
                    sg=0.001*ones(N,1);
                    for t=0:dt:Tmax

                        % Adding in a lot of noise to V1 by itself

                        if(seed>0)
                            Istim(seed) = 0.365+0.005*randn;
                            Istim(seed+N/2) = 0.365+0.005*randn;
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
%                     keyboard
                    for region=1:N
                        data=neuro_act2(sub,1:end,region)';
                        [acf,lags] = autocorr(data(5e3:end),5000);
                        acf_data(:,region,sub) = acf;
                    end
                    
                end
                
            end

        end
        acf=mean(acf_data,3);
        save(filename,'acf','SEED','sigma');
end

