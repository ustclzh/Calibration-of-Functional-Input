function RES_opt=simulation_for_one_basis(PRE,True,prior_specification,mcmc)
% PRE： contains the information for basis decomposition, initial design and emulator;
% True: information for true value of functional input and observation data;
% prior_specification: indicate the prior model used in simulation, 1: Unif-GP; 2: L-GP; 3: M-GP; 4: U-GP;
% design_criterion: indicate the design criterion used in simulation, 1: WPV; 2: VL; 3: AEI;
% design_size: the size of follow-up design, if set as 0, means use initial design directly;
% mcmc: can be empty, if not empty, the simulation will run mcmc algorithm.
%%
pre=PRE{prior_specification};
data=pre.data;
data.obs=True.obs;
data.sigma_e=True.sigma_e;
%%
start_set=randn(20,size(data.KL,1));
result=H2P_simulator(data,start_set);
z_map=result.x_opt(1:data.M);
true=True.f_input;
result.post_err=mean((z_map*data.KL'-true).^2)^(1/2)/(mean(true.^2)^(1/2));
RES_opt{1}=result;
%% MCMC part
if nargin==4
start_point=[result.x_opt];
para_MCMC.nsimu=60000;
para_MCMC.burnin=10001;
para_MCMC.fixeta=data.fix_eta;% 0 for random eta, 1 for fixed eta
if para_MCMC.fixeta==1
    para_MCMC.eta=data.eta;
end
para_MCMC.sigmasq_f_upper=15;
results=MwG_Simulator(data,True,para_MCMC,start_point);
results.post_err=mean((results.f_postmean'-results.True.f_input).^2)^(1/2)/mean((results.True.f_input).^2)^(1/2);
RES_opt{2}=results;
end
end
%%