clear,clc
warning off;
global nugterm
nugterm=10^(-6);
addpath(genpath(pwd));
run('prepare.m');
sigma_e=0.02;
data=pre_KL.data;
PRE={pre_KL,pre_full};
eta=[0.01*ones(1,30),0.5*ones(1,30),0.99*ones(1,30)];
%%
for i=1:90
    i
    eta_true=eta(i);%rand(1)*0.98+0.01;
    Z_true=randn(1,21);
    [True,~]=Pde_parameter(eta_true,Z_true,sigma_e,1);
    design_size=20;
    design_criterion=1;
    RES_opt{1}=simulation_for_one_basis_one_design(PRE,True,1,design_criterion,design_size);
    ESS_post(i,1)=RES_opt{1}{1}.post_err;
    RES_opt{2}=simulation_for_one_basis(PRE,True,1);
    ESS_post(i,2)=RES_opt{2}{1}.post_err;
    RES_opt{3}=simulation_for_one_basis(PRE,True,2);
    ESS_post(i,3)=RES_opt{3}{1}.post_err;
    RES_OPT{i}=RES_opt;
end
filename=['RESULT_exact.mat'];
save(filename,'RES_OPT','-v7.3')
filename=['ESS_exact.mat'];
save(filename,'ESS_post','-v7.3')


