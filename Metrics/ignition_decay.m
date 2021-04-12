function [Ignition2,Excitability2,ignition1,ignition1a,recalcuated_decay] = ignition_decay(PERTURB,neuro_act_all)   
    % Parameters that were used in the simulation. 
    kk=1;
    N=68;
    tseed=1:N;
    SEED=10;

    % First assign the variables to get the max peak rate:
    for perturb=PERTURB
    	neuro_actf=neuro_act_all(:,:,kk);
    	ssnum=1;
	    for ss=tseed
	        peakrate3(ssnum)=mean(neuro_actf(160:200,ss),1);
	        basal(ssnum)=mean(neuro_actf(100:140,ss));
	        tscale=0:1/length(decayneuro):1;
	        ssnum=ssnum+1;	        
	    end
        peakrate2(kk,:)=peakrate3./basal;
	    kk=kk+1;
	end
	
	% Below we have the range that you want the function to eventually sample to:
	PERTURB_sample=0.:0.001:0.2;    

	% Fitting rourtines for ignition
    sigfunc = @(A, x)(A(1) ./ (1 + exp(-A(2)*(x-A(3)))) + A(4));
	options=optimset('MaxFunEvals',10000,'MaxIter',1000,'Display','off');
    
    peakrate1=mean(peakrate2(end-10,:));
    nntarget=1;
    figure('color','white');
    for ntarget=1:N
        A0=[mean(peakrate2(end-10:end,ntarget))-mean(peakrate2(1:10,ntarget)) 10 0.1 mean(peakrate2(1:10,ntarget))];
        Afit = lsqcurvefit(sigfunc,A0,PERTURB,(peakrate2(:,ntarget))',[0 0 -1 0],[100 100 1 10*mean(peakrate2(1:10,ntarget))],options);        
        % Just a little part here to sample the fitted regime to a different level (the normalization level is always for above)
        yfit=Afit(1) ./ (1 + exp(-Afit(2)*(PERTURB_sample-Afit(3))))+Afit(4);
        if(ntarget<35)
        	subplot(6,6,ntarget)
        	plot(PERTURB,(peakrate2(:,ntarget))','.');hold on;plot(PERTURB_sample,yfit)
        	title(['Region: ',num2str(ntarget)])
        	xlabel('$\rho$','Interpreter','LaTeX');ylabel('$r^E_{max}/ r^E_{rest}$','Interpreter','LaTeX');
        	legend({'Data','fit'})
        end
        % Calculating ignition and remembering to divide by the sampling rate (it shoudl really just be "grad" to make it easier..)
        ignition1(nntarget)=max(diff(diff(yfit/0.001)/0.001))*yfit(end)/1000;
        ignition1a(nntarget)=yfit(end)/1000;
        ignition1b(nntarget)=max(diff(diff(yfit/0.001)/0.001));
        nntarget=nntarget+1;
    end

    % Here we are looking at the SEED 10, and its contralateral component and zeroing it out
    ignition1(SEED)=0;
    ignition1(SEED+N/2)=0;

    Ignition2=mean(ignition1(find(ignition1>mean(ignition1)+std(ignition1))));
    Excitability2=mean(peakrate1);    


    % Now here we want to calculate the decay using a nonlinear fit 
    kk=size(neuro_act_all,3);
    % Using the last point (the point of greatest stimulation)
    sr = 20*1e-3;
	time_vector=120:350;
	time = 0:sr:(length(time_vector)-1)*sr;
	whole_time=0:sr:(350-1)*sr;

	% Using a nonlinear decay function instead
	s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[1 0.5 3]);
	f = fittype('A*(exp(-x*D)+B)','options',s);
	new_time=linspace(whole_time(1),whole_time(end),1000);
	stim_vector = perturb*0.5*(sign(new_time-3) - sign(new_time-4));



	neuro_actf=neuro_act_all(:,:,kk);
	ssnum=1;
    for ss=tseed
        decayneuro=squeeze(neuro_actf(201:end,ss));
    	tscale=0:sr:(length(decayneuro)-1)*sr;
    	[c, gof] = fit(tscale(1:length(decayneuro))',decayneuro,f);
    	recalcuated_decay(ss)=c.D;  
    	ssnum=ssnum+1;
    	% figure;    
	    % plot(tscale(1:length(decayneuro)),(decayneuro'))
	    % hold on;
	    % plot(tscale(1:length(decayneuro)),f(c.A,c.B,c.D,tscale(1:length(decayneuro))))
	    % plot(tscale(1:length(decayneuro)),exp(polyval(bdecay,tscale(1:length(decayneuro)))));
	    % legend({'\phi(t)','fitNew','fitOLD'});
    end
end
