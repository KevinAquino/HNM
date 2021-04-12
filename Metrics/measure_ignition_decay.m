% This function loads all up the pre-calculated simulations and calculates metrics from them

modelType = {'EI','PCA','T1T2','BEI','Null'};
simulations=[1:999];
Perturb_range=[1];

% First load up the metric

for model=1:length(modelType),
    duds=[];
	for sim=simulations,
        disp(['UP TO......',modelType{model},' sim ',num2str(sim),'/999']);
        try
            for pert = Perturb_range,
                fileName=['/scratch/kg98/kevo/Ignition_decay/',modelType{model},'/',modelType{model},'_',num2str(sim),'_Peturb_range_',num2str(pert),'.mat'];
                % Load up the file
                load(fileName);
                % Now calculate the metric - do some agglomeration if need be
            end
            [Ignition2,Excitability2,ignition1,ignition1a,recalcuated_decay] = ignition_decay(PERTURB,neuro_act_all);
        catch
            ignition1=0*ignition1;ignition1a=0*ignition1a;
            duds = [duds,sim];
        end
		% Saving all the responses:
		ignition_all(model,sim,:)=ignition1;
		response_all(model,sim,:)=ignition1a;
		ignition_var_all(model,sim,:)=mean(ignition1(find(ignition1>mean(ignition1)+std(ignition1))));        
    end
    remove_list{model}=duds;
    save all_data
end