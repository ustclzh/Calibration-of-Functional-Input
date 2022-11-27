clear,clc
warning off;
addpath(genpath(pwd));
run('prepare.m');
%load('pre9basis.mat');
%load('inputdesign.mat');
PRE={pre_KL_aniso,pre_KLL,pre_KLC,pre_KLU};
data=pre_KL_aniso.data;
eta=[data.eta_L*ones(1,10),0.5*(data.eta_L+data.eta_U)*ones(1,10),data.eta_U*ones(1,10)];
K=90;
bar = waitbar(0,'Initializing');
for i=1:K
str=['Compare Criteria',num2str(i),' of ' num2str(K)];   
waitbar(i/K,bar,str)  
eta_true=eta(rem(i-1,30)+1);
Z_true=design_z(i,:);
sigma_e=0.02;
[True,~]=Pde_parameter(eta_true,Z_true,sigma_e,1);
for design_criterion=1:3
    basis_type=1;
    design_size=20;
    RES_opt{design_criterion}=simulation_for_one_basis_one_design(PRE,True,basis_type,design_criterion,design_size);
    ESS_post(i,design_criterion)=RES_opt{design_criterion}{1}.post_err;
end
RES_OPT{i}=RES_opt;
end
%
filename=['RESULT_criteria','.mat'];
save(filename,'RES_OPT','-v7.3')
filename=['ESS_criteria','.mat'];
save(filename,'ESS_post','-v7.3')

