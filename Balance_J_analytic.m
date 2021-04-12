function [Joptim]=Balance_J_analytic(we,C,gain)

	gain=gain(1,:)';
	dtt   = 1e-3;   % Sampling rate of simulated neuronal activity (seconds)
	dt=0.1;
	taon=100;
	taog=10;
	gamma=0.641;	
	JN=0.15;
	I0=0.382;
	Jexte=1.;
	Jexti=0.7;
	w=1.4;
	r_0=3.1;

	% keyboard
	% A rescaling here to make it easier.
	taon=taon/1000;
	taog=taog/1000;


	% Solved by working out the transcendtal equation:
	ii=linspace(0,1,10000);	
	% keyboard
	functional=Jexti*I0 + JN*gamma*r_0*taon/(1+ gamma*r_0*taon) - taog*phii_gain(ii,gain) - ii;
	[~,ind]=min(abs(functional));
	ii=ii(ind);

	disp(['Solved transcendtal equation numerically, Ii(I)=',num2str(ii)]);

% keyboard
	Joptim = sum(C)'*we*JN/(phii_gain(ii,gain)*taog)*(gamma*r_0*taon/(1+ taon*gamma*r_0)) + 1;

	% Transpose for convention
	% Joptim=Joptim';
