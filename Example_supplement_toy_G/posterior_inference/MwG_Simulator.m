function [result,chain]=MwG_Simulator(data,True,para_MCMC,start_point)
global gradeiant_G chain
% para_emulator: all parameters for the emulator trained bt a design,
% para_MCMC.nsimu: length of MCMC
% start_point: starting point of the chain
% para_MCMC.fixeta: an indicator for whether fixing eta. 1 for fixed 0 for not fixed
% para_MCMC.eta : the eta value that is fixed. Only valid when fixeta=1.
%% preparation
scale=data.M/2.4^2;
time_used=0;
n_simu=para_MCMC.nsimu;
M=data.M;
data.sigma_e=True.sigma_e;
Z=start_point(1:M);
if para_MCMC.fixeta==1  % assume eta known
    p=M+1;
   eta=data.eta;
   result.fixeta=eta;
   invR_Z_eta=eye(M);
   R_Z_eta=eye(M);
else % assume eta random
    p=M+2-data.isotropic;
    eta=0.5*(data.eta_L+data.eta_U)*ones(1,2-data.isotropic);
    invR_Z_Eta=data.invR_Z_Eta;
    logdetR_Z_Eta=data.logdetR_Z_Eta;
    R_Z_Eta=data.R_Z_Eta;
end
emulator_Z=Simulator(Z,data);
Grad_Z=emulator_grad(Z,data);%gradientofemulator(Z,emulator_Z,data,para_emulator);
chain=zeros(n_simu,p+1);
Lik=zeros(1,n_simu);

Err1=zeros(1,n_simu);
Err2=zeros(1,n_simu);
Err_out=zeros(1,n_simu);
Errinfty=zeros(1,n_simu);
h=waitbar(0,'start');
accepted_num=0;
sigmasqf=1;
%start the MCMC
%mu=start_point(1:M)*True.KL';
for i=1:n_simu
    tic
    if para_MCMC.fixeta==0    % if eta random, update eta;
         if data.isotropic==1
            eta=slicesampling(eta,Z,invR_Z_Eta,logdetR_Z_Eta,data,sigmasqf);
            ind=floor((eta-data.eta_L)*10000+1.0001);
            invR_Z_eta=invR_Z_Eta{ind};
            R_Z_eta=R_Z_Eta{ind};
        else
            eta=slicesampling_multi(eta,Z,invR_Z_Eta,logdetR_Z_Eta,data,sigmasqf);
            ind1=floor((eta(1)-data.eta_L)*100+1.0001);
            ind2=floor((eta(2)-data.eta_L)*100+1.0001);
            invR_Z_eta=invR_Z_Eta{ind1,ind2};
        end
    end
    [Z_new,sigma_z_current]=proposal(Z,invR_Z_eta,Grad_Z,data,sigmasqf,R_Z_eta,scale); % a new proposal
    emulator_Z_new=Simulator(Z_new,data); %
    gradeiant_G=Grad_Z;
    if rem(i,1000)==1
        Grad_Z_new=emulator_grad(Z_new,data);%gradientofemulator(Z_new,emulator_Z_new,data,para_emulator);
    else
        Grad_Z_new=Grad_Z;
    end
    % reject or accept
    logpold2new=logtransform(Z_new,Z,invR_Z_eta,Grad_Z,data,sigmasqf,R_Z_eta,scale);
    logpnew2old=logtransform(Z,Z_new,invR_Z_eta,Grad_Z_new,data,sigmasqf,R_Z_eta,scale);% these two are log transform density of old to new and new to old.
    logpznew=logconditionalZ(Z_new,emulator_Z_new,invR_Z_eta,data,sigmasqf,R_Z_eta); % log posterior of new proposal
    Lik(i)=logpznew;
    logpz=logconditionalZ(Z,emulator_Z,invR_Z_eta,data,sigmasqf,R_Z_eta);% log posterior of current state.
    Pr=exp((logpznew-logpz+logpnew2old-logpold2new));
    u=rand(1);
    if Pr>=u %accept new proposal
        Z=Z_new;
        accepted_num=accepted_num+1;
        Grad_Z=Grad_Z_new;
        emulator_Z=emulator_Z_new;
    end
    sigmasqf=conditionalsigmasq_f(invR_Z_eta,Z_new,data);
    chain(i,:)=[Z,eta,sigmasqf]; % new sample
    Err1(i)=(mean(abs(Z*data.KL'-True.f_input)));
    Err_out(i)=(mean((emulator_Z-True.obs).^2))^(1/2);
    Err2(i)=(mean((Z*data.KL'-True.f_input).^2)/mean((True.f_input).^2))^(1/2);
    Errinfty(i)=(max(abs(Z*data.KL'-True.f_input)));
    time_Z=toc;
    time_used=time_used+time_Z;
    if rem(i,500)==0
    str=[num2str(floor(n_simu/1000-i/1000)) 'K steps, ' num2str(floor(time_used*(n_simu-i)/i)), 's.',' P_a: ' num2str(floor(100*accepted_num/i))];
    waitbar(i/n_simu,h,str);
    end
end
close(h);
result.likelihood=Lik;
result.data=data;
result.True=True;
result.Err_out=Err_out;
result.Err2=Err2;
result.Err1=Err1;
result.Errinfty=Errinfty;
result.Z_MAP=MAP_est(chain((n_simu/6+1):50:end,1:data.M));% MAP estimation of KL coefficients
%result.obs=Obs;
result.chain=chain;%MCMC result
result.accepted_rate=accepted_num/n_simu; % acceptance rate of the Metropolis-Hastings algorithm for sampling from the full conditional distribution of \zeta
jump=floor((para_MCMC.nsimu-para_MCMC.burnin+1)/2500);
result.Z_postmean=mean(chain(para_MCMC.burnin:jump:end,1:data.M));% posterior mean of KL coefficients
result.f_postmean=data.KL*result.Z_postmean';% posterior mean of functional input
result.time=sum(time_used);% cpu time used for the whole MCMC
result.para_MCMC=para_MCMC;
disp(['MCMC algorithm time ' num2str(time_used) 'seconds']);
end
%% functions
%% slice sampling
function y=slicesampling(eta_old,Z,invR,logdetR,data,sigmasqf) 
eta_L=data.eta_L;
eta_U=data.eta_U;
% slice sampling 
for i=1:10
    u=log(rand(1))+logconditionaleta(eta_old,Z,invR,logdetR,data,sigmasqf);
    L=eta_L;
    R=eta_U;
    temp=rand(1)*(R-L)+L;
    u_temp=logconditionaleta(temp,Z,invR,logdetR,data,sigmasqf);
    while u_temp<u
        if temp<eta_old
           L=max(L,temp); 
        end
        if temp>eta_old
           R=min(R,temp);            
        end
        temp=rand(1)*(R-L)+L;
        u_temp=logconditionaleta(temp,Z,invR,logdetR,data,sigmasqf);
    end
    if temp<eta_L
        eta_new=eta_L;
    elseif temp>eta_U
        eta_new=eta_U;
    else
        eta_new=temp;
    end
    eta_old=eta_new;
end
y=eta_new;
end
function y=logconditionaleta(eta,Z,invR,logdetR,data,sigmasqf)% using slise sampling
ind=round((eta-data.eta_L)*9800/(data.eta_U-data.eta_L)+1);
invSigma_z=invR{ind}/sigmasqf;
y=-Z*(invSigma_z*Z')/2-logdetR(ind)/2;
end
%%
function y=slicesampling_multi(eta_old,Z,invR,logdetR,data,sigmasqf) 
eta_L=data.eta_L;
eta_U=data.eta_U;
% slice sampling 
for j=1:2
    eta_temp=eta_old;
   for i=1:10
    u=rand(1)*conditionaleta_multi(eta_old,Z,invR,logdetR,data,sigmasqf);
    L=eta_L;
    R=eta_U;
    temp=rand(1)*(R-L)+L;
    eta_temp(j)=temp;
    u_temp=conditionaleta_multi(eta_temp,Z,invR,logdetR,data,sigmasqf);
    while u_temp<u
        if temp<eta_old(j)
           L=max(L,temp); 
        end
        if temp>eta_old(j)
           R=min(R,temp);            
        end
        temp=rand(1)*(R-L)+L;
        eta_temp(j)=temp;
        u_temp=conditionaleta_multi(eta_temp,Z,invR,logdetR,data,sigmasqf);
    end
    eta_new=eta_temp;
    eta_old=eta_new;
   end
    if eta_new(j)<eta_L
        eta_new(j)=eta_L;
    elseif eta_new(j)>eta_U
        eta_new(j)=eta_U;
    end
    y=eta_new;
end
end
function y=conditionaleta_multi(eta,Z,invR,logdetR,data,sigmasqf)
ind1=round((eta(1)-data.eta_L)*100+1);
ind2=round((eta(2)-data.eta_L)*100+1);
invSigma_z=invR{ind1,ind2}/sigmasqf;
y=exp(-Z*(invSigma_z*Z')/2)*(exp(logdetR(ind1,ind2)))^(-1/2); 
end
%%
function y=conditionalsigmasq_f(invR_zeta_eta, Z,data) 
%return an inverse gamma sample
% the prior is invgamma(a,b)
a=data.para_prior.a;
b=data.para_prior.b;
a=a+data.M/2;
b=b+Z*(invR_zeta_eta*Z')/2;
y=(gamrnd(a,1/b))^(-1);
if y>100
    y=100;
end
if isnan(y)
    y=100;
end
end
%%
function y=logconditionalZ(Z,emulator_Z,invR_Z_eta,data,sigmasqf,R_Z_eta)
Sigma_e=data.sigma_e^2;
logprior=-Z*(R_Z_eta\Z')/(2*sigmasqf);
loglikelihood=-sum((data.obs-emulator_Z).^2)/(2*Sigma_e);
y=loglikelihood+logprior; %logposterior
end
function [y,sigma]=proposal(Z,invR_Z_eta,Grad,data,sigmasqf,R_Z_eta,scale) % generating new proposal
M=data.M;
Sigma_e=data.sigma_e^2;
sigma=invandlogdet(Grad*Grad'/Sigma_e+invR_Z_eta/sigmasqf);
sigma=sigma/scale;
y=mvnrnd(zeros(1,data.M),sigma)+Z;
end
function y=logtransform(Z_new,Z,invR_Z_eta,Grad,data,sigmasqf,R_Z_eta,scale) 
M=data.M;
Sigma_e=data.sigma_e^2;
[sigma,logdetinvsig]=invandlogdet(Grad*Grad'/Sigma_e+invR_Z_eta/sigmasqf);
sigma=sigmasqf*R_Z_eta-sigmasqf*R_Z_eta*Grad*(data.sigma_e^2*eye(size(data.obs,2))+Grad'*R_Z_eta*Grad)*Grad'*R_Z_eta;
sigma=sigma/scale;
y=-(Z_new-Z)*((sigma)\(Z_new'-Z'))/2+logdetinvsig/2;
end
function y_grad=emulator_grad(x,data)
n2=length(data.obs);
origin=Simulator(x,data); 
R3=zeros(length(x),n2);
for j=1:length(x)
    x_p=x;
    x_p(j)=x_p(j)+0.1;
    R3(j,:)=(Simulator(x_p,data)-origin)/0.1;
end
y_grad=R3;
end
function y=MAP_est(sample)
    %return the MAP estimation for parameter, each column stands for a group of sample for one parameter, return the density mode of each
    %column.
    n=size(sample,2);
    for i=1:n
    data=sample(:,i);
    [a,b]=ksdensity(data);
    a_max=find(a==max(a));
    y(i)=b(a_max(1));
    end
end