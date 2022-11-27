function [True,data]=Pde_parameter(eta_true,z_true,sigma_e,isotropic)
% check
% something need to be changed when pde model change
% this function is for generating the true value, observations, determine
% the PDE geometry and parameters.
%sigma_e the propotion of error added into the observations.
% all global variables are PDE parameters.
data.eta_L=0.01;
data.eta_U=0.99;
if nargin>3
    data.isotropic=isotropic;%0 means input GP is anisotropic, 1 means input GP is isotropic.
else
    data.isotropic=1;
end
X=(0:20)'/20;
Y=zeros(size(X,1),1);
True.X=X;
True.Y=Y;
True.z_true=z_true;
data.D_in=(1:1:19)'/20;
data.t_input=0:0.02:0.2;
True.t_input=0:0.02:0.2;
%% Parameter setting
True.eta_true=eta_true;
R=@(X,Y,eta) exp(log(eta)*(((X(:)-X(:)').^2+(Y(:)-Y(:)').^2)))+eye(size(X,1))*10^(-6);
True.z_true=z_true;
True.sigma_e=sigma_e;
True.R=R;
True=PRE_SETUP(True);
%data.R_int=R_int;
data.obs=True.obs;
data.X=True.X;
data.Y=True.Y;
data.sigma_e=sigma_e;
data.R=R;
data.explain_var=cumsum(True.eig)/sum(True.eig);
end

function True=PRE_SETUP(parameter_setting)
True=parameter_setting;
Sigma=True.R(True.X,True.Y,True.eta_true);
[S,V,~]=svd(Sigma);
M_true=find(cumsum(diag(V))./sum(diag(V))>0.95, 1);
N=length(True.X);
KL_true=S*(V.^0.5);%each column is the principal components. 
z_true=True.z_true;
M=length(z_true);
True.KL=KL_true;
if N>M
    z_true(M+1:N)=zeros(1,N-M);
end
if N<M
    z_true=z_true(1:N);
end
f_input=z_true*KL_true';
True.eig=diag(V);
True.M=M_true;
True.f_input=f_input;
True.z_true=z_true;
True.obs_no_err=Simulator(z_true,True);
True.obs=True.obs_no_err+True.sigma_e*randn(1,length(True.obs_no_err));
end