clear,clc
warning off;
addpath(genpath(pwd));
run('prepare391215.m');
%load('pre391215basis.mat');
%load('inputdesign.mat');
PRE={pre_KL3,pre_KL9,pre_KL12,pre_KL15};
data=pre_KL9.data;
eta=[data.eta_L*ones(1,10),0.5*(data.eta_L+data.eta_U)*ones(1,10),data.eta_U*ones(1,10)];
K=30;
bar = waitbar(0,'Initializing');
for i=1:K
str=['Compare Criteria',num2str(i),' of ' num2str(K)];   
waitbar(i/K,bar,str)  
eta_true=eta(i);
Z_true=design_z(i,1:10);
sigma_e=0.02;
[True,~]=Pde_parameter(eta_true,Z_true,sigma_e,1);
for basis_type=1:4
    design_size=20;
    design_criterion=1;
    RES_opt{basis_type}=simulation_for_one_basis_one_design(PRE,True,basis_type,design_criterion,design_size);
    ESS_post(i,basis_type)=RES_opt{basis_type}{1}.post_err;
    ESS_prior(i,basis_type)=RES_opt{basis_type}{1}.prior_err;
end
RES_OPT{i}=RES_opt;
end
close(bar)
filename=['RESULT_M','.mat'];
save(filename,'RES_OPT','-v7.3')
filename=['ESS_M','.mat'];
save(filename,'ESS_post','-v7.3')

