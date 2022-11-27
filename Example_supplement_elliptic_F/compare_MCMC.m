clear,clc
warning off;
addpath(genpath(pwd));
run('prepare.m');
%load('pre9basis.mat');
%load('inputdesign.mat');
PRE={pre_KL,pre_KLL,pre_KLC,pre_KLU};
data=pre_KL.data;
eta=[data.eta_L*ones(1,10),0.5*(data.eta_L+data.eta_U)*ones(1,10),data.eta_U*ones(1,10)];
eta_true=0.5;
Z_true=design_z(1,:);
sigma_e=0.02;
[True,~]=Pde_parameter(eta_true,Z_true,sigma_e,1);
for basis_type=1:4
    design_size=20;
    design_criterion=1;
    RES_opt{basis_type}=simulation_for_one_basis_one_design(PRE,True,basis_type,design_criterion,design_size,1);
    ESS_post(basis_type)=RES_opt{basis_type}{1}.post_err;
end
filename=['RESULT_MCMC','.mat'];
save(filename,'RES_opt','-v7.3')
filename=['ESS_MCMC','.mat'];
save(filename,'ESS_post','-v7.3')