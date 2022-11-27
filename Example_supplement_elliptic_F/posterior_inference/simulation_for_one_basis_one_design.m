function RES_opt=simulation_for_one_basis_one_design(PRE,True,prior_specification,design_criterion,design_size,mcmc)
% PRE： contains the information for basis decomposition, initial design and emulator;
% True: information for true value of functional input and observation data;
% prior_specification: indicate the prior model used in simulation, 1: Unif-GP; 2: L-GP; 3: M-GP; 4: U-GP;
% design_criterion: indicate the design criterion used in simulation, 1: WPV; 2: VL; 3: AEI;
% design_size: the size of follow-up design, if set as 0, means use initial design directly;
% mcmc: can be empty, if not empty, the simulation will run mcmc algorithm.
%%
pre=PRE{prior_specification};
data=pre.data;

data.Candidate_for_start=pre.Candidate_for_start;
para_emulator=pre.para_emulator;
Candidate_for_start=pre.Candidate_for_start;
data.obs=True.obs;
data.obsp=data.obs*para_emulator.s;
data.sigma_e=True.sigma_e;
n_batch_start=data.M*3;
start_set=generate_design_support_points(n_batch_start,Candidate_for_start,100);
%%
result=H2P(data,para_emulator,start_set);
z_initial=result.x_opt(1:data.M);
para_design.initial_design=para_emulator.D_para;
para_design.cor_func=para_emulator.cor_func;
para_design.D_in=para_emulator.D_in;
para_design.d_desire=design_size; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% design size
[para_design_adv,para_emulator_adv]=fully_sequential_generation(data,para_emulator,para_design,z_initial,design_criterion);
Candidate_for_start=[start_set;para_design_adv.D_para];
start_set=generate_design_support_points(n_batch_start,Candidate_for_start,100);
start_set=[start_set;z_initial];
[para_emulator_adv,data]=Train_Emulator(data,para_design_adv,para_emulator_adv);
result=H2P(data,para_emulator_adv,start_set);
f_true=True.f_input;
z_project=((data.KL'*data.KL)\(data.KL'*f_true'))';
z_map=result.x_opt;
f_map=z_map(1:data.M)*data.KL';
relativel2e=(mean((f_true-f_map).^2))^(1/2)/(mean(f_true.^2)^(1/2));
[m,~]=emulator(z_map(1:data.M),para_emulator_adv);
RSS=mean((m-data.obs).^2)^(1/2);
[m,~]=emulator(z_project,para_emulator_adv);
RSS_o=mean((m-data.obs).^2)^(1/2);
result.post_err=relativel2e;
result.z_project=z_project;
result.z_true=True.z_true;
result.prior_err=mean((z_project*data.KL'-f_true).^2)^(1/2)/(mean(f_true.^2)^(1/2));
result.f_post=f_map;
result.para_design=para_design_adv;
result.obs=data.obs;
result.eta_true=True.eta_true;
%result.True=True;
result.f_true=f_true;
result.f_proj=z_project*data.KL';
RES_opt{1}=result;
[RSS,RSS_o,result.prior_err,relativel2e];%list result
%% MCMC part
if nargin==6
start_point=result.x_opt;
para_MCMC.nsimu=60000;
para_MCMC.burnin=10001;
para_MCMC.fixeta=data.fix_eta;% 0 for random eta, 1 for fixed eta
if para_MCMC.fixeta==1
    para_MCMC.eta=data.eta;
end
para_MCMC.multi_design=design_size;
para_MCMC.sigmasq_f_upper=15;
restrict_emu=10;
results=MwG(data,True,para_emulator_adv,para_MCMC,restrict_emu,start_point);
RES_opt{2}=results;
end
end
%%
function [y_pred,PosVar]=emulator(x,para_emulator)
D_in=para_emulator.D_in;
D_para=para_emulator.D_para;
cor_func=para_emulator.cor_func;
thetaopt=para_emulator.thetaopt;
alpha=para_emulator.alpha;
sigma2=para_emulator.sigma2;
%R1=para_emulator.R1;
ResR2I=para_emulator.ResR2I;
R2I=para_emulator.R2I;
%MTRIMI=para_emulator.MTRIMI;
M2TRI=para_emulator.M2TRI;
[n1, d1]=size(D_in);
[n2, d2]=size(D_para);
% R3=zeros(1,n2);
% for j=1:n2
%     R3(j)=CompCorr(x,D_para(j,:),thetaopt((d1+1):end),cor_func);
% end
R3=CompCorr_1(x,D_para,thetaopt((d1+1):end));
Ypred=alpha+ResR2I*(R3');
y_pred=Ypred';
PosVar=max(sigma2*(1-R3*R2I*R3')+sigma2*(1-M2TRI*R3')^2/(sum(M2TRI)),0);
end
function y=CompCorr_1(x1,D,eta)
    n=size(D,1);
    rho1=sqrt(6)*abs(x1-D)./(ones(n,1)*eta);
    y=prod((exp(-rho1')).*(rho1'+1));
end

%%
function y=generate_design_support_points(n,Sample,maxIte)
% warning: chack that Sample doesn't contain any repeatation
% support points
% n: desired design size
% Sample: candidate sample
% maxIte: max iteration, default 200;

D_new=Sample(1:n,:);
Sample(1:n,:)=[];
p=size(Sample,2);
D_0=zeros(n,p);
if nargin==2
    maxIte=200;
end
ite=0;
while sum(abs(D_0-D_new))>0.01
    ite=ite+1;
    D_0=D_new;
   for i=1:n
       D_new(i,:)=Mi(D_0,i,Sample);
   end
   if ite>=maxIte
       break;
   end
end
y=D_new;
end
function y=Mi(D,i,sample)
x=D(i,:);
N=size(sample,1);
p=length(x);
n=size(D,1);
d_2=sum((ones(N,1)*x-sample).^2,2).^(1/2);
q=sum(d_2.^(-1));
D(i,:)=[];
d_x=sum((ones(n-1,1)*x-D).^2,2).^(1/2);
D_w=(x-D)./(d_x*ones(1,p));
x_new1=N*sum(D_w,1)/n;
x_new2=sum(sample./(d_2*ones(1,p)),1);
y=(x_new1+x_new2)/q;
end
