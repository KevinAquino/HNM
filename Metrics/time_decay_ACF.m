function taus = time_decay_ACF(acf)

% Fitting functions:
s_2 = fitoptions('Method','NonlinearLeastSquares','StartPoint',[0.3 0.6 0.1 0.2],'Lower',[0 0 0 0]);
f_2 = fittype('A*exp(-x*D) + B*exp(-x*D2)','options',s_2);


s_1 = fitoptions('Method','NonlinearLeastSquares','StartPoint',[1 0.6],'Lower',[0 0]);
f_1 = fittype('A*exp(-x*D)','options',s_1);


% Below here are the plots for fitting
for region=1:68,
		acinds=find(acf(:,region)>0.05);
		% Assume ACF is monotonically decreasing
		acf_use=acf(1:max(acinds)',region);
		time_use=1e-3*(0:(max(acinds)-1));
		% time_use=20e-3*(0:(max(acinds)-1));
		if(acinds<10)
			taus(:,region)=1e-3;
		else
			% [c, gof] = fit(time_use',acf_use,f);
			[c1, gof1] = fit(time_use',acf_use,f_1);
			[c2, gof2] = fit(time_use',acf_use,f_2);			
			% Here testing which model is better, if sse of single is 8 times less than that of the double report it other wise do the double.
			if(gof1.sse/gof2.sse<8)
				c=c1;
				taus(1,region)=1./c.D;
				taus(2,region)=1./c.D;
				taus(3,region)=1./c.D;		
			else
				c=c2;
				taus(1,region)=1./c.D;
				taus(2,region)=1./c.D2;
				taus(3,region)=1/2*(c.A./c.D + c.B./c.D2);							
		end
	end
end