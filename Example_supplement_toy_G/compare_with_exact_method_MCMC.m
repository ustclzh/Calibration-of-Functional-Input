clear,clc
warning off;
global nugterm
nugterm=10^(-6);
addpath(genpath(pwd));
run('prepare.m');
sigma_e=0.02;
data=pre_KL.data;
PRE={pre_KL,pre_full};
eta_true=0.4;
Z_true=randn(1,21);
[True,~]=Pde_parameter(eta_true,Z_true,sigma_e,1);
design_size=20;
design_criterion=1;
RES_opt{1}=simulation_for_one_basis_one_design(PRE,True,1,design_criterion,design_size,1);
RES_opt{2}=simulation_for_one_basis(PRE,True,1,1);
RES_opt{3}=simulation_for_one_basis(PRE,True,2,1);
for j=1:3
    ESS_post(j)=RES_opt{j}{1}.post_err;
end
filename=['RESULT_MCMC.mat'];
save(filename,'RES_opt','-v7.3')
filename=['ESS_MCMC.mat'];
save(filename,'ESS_post','-v7.3');

