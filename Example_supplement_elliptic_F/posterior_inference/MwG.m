function [result,chain]=MwG(data,True,para_emulator,para_MCMC,restrict_emu,start_point)
% para_emulator: all parameters for the emulator
% para_MCMC.nsimu: length of MCMC
% start_point: starting point of the chain
% para_MCMC.fixeta: an indicator for whether fixing eta. 1 for fixed 0 for not fixed
% para_MCMC.eta : the eta value that is fixed. Only valid when fixeta=1.
%% preparation
global scale
scale=data.M/2.4^2;
time_used=0;
n_simu=para_MCMC.nsimu;
M=data.M;

Z=start_point(1:M);
if para_MCMC.fixeta==1  % assume eta known
    p=M+1;
   eta=data.eta;
   result.fixeta=eta;
   invR_Z_eta=eye(M);
   
else % assume eta random
    p=M+2-data.isotropic;
    eta=0.5*ones(1,2-data.isotropic);
    invR_Z_Eta=data.invR_Z_Eta;
    logdetR_Z_Eta=data.logdetR_Z_Eta;
end
[emulator_Z,Cov_post]=emulator_trans(Z,para_emulator);
Grad_Z=emulator_grad(Z,para_emulator);%gradientofemulator(Z,emulator_Z,data,para_emulator);
chain=zeros(n_simu,p+1);
Lik=zeros(1,n_simu);
Var_pred=zeros(n_simu,1);
Err1=zeros(1,n_simu);
Err2=zeros(1,n_simu);
Errinfty=zeros(1,n_simu);
h=waitbar(0,'start');
accepted_num=0;
[S_R1,V_R1,D_R1]=svd(para_emulator.R1);
data.obsp=data.obs*S_R1;
sigmasqf=1;
%start the MCMC
for i=1:n_simu
    tic
    if para_MCMC.fixeta==0    % if eta random, update eta;
         if data.isotropic==1
            eta=slicesampling(eta,Z,invR_Z_Eta,logdetR_Z_Eta,data,sigmasqf);
            ind=floor((eta-data.eta_L)*10000+1.0001);
            invR_Z_eta=invR_Z_Eta{ind};
        else
            eta=slicesampling_multi(eta,Z,invR_Z_Eta,logdetR_Z_Eta,data,sigmasqf);
            ind1=floor((eta(1)-data.eta_L)*100+1.0001);
            ind2=floor((eta(2)-data.eta_L)*100+1.0001);
            invR_Z_eta=invR_Z_Eta{ind1,ind2};
 
        end
    end
    [Z_new,sigma_z_current]=proposal(Z,invR_Z_eta,Grad_Z,data,Cov_post,S_R1,V_R1,sigmasqf); % a new proposal
    [emulator_Z_new,Cov_post_new]=emulator_trans(Z_new,para_emulator); %
    Var_pred(i)=Cov_post_new;
    %Output(i,:)=emulator_Z_new;
    Grad_Z_new=emulator_grad(Z_new,para_emulator);%gradientofemulator(Z_new,emulator_Z_new,data,para_emulator);
    % reject or accept
    logpold2new=logtransform(Z_new,Z,invR_Z_eta,Grad_Z,data,Cov_post,S_R1,V_R1,sigmasqf);
    logpnew2old=logtransform(Z,Z_new,invR_Z_eta,Grad_Z_new,data,Cov_post_new,S_R1,V_R1,sigmasqf);% these two are log transform density of old to new and new to old.
    logpznew=logconditionalZ(Z_new,emulator_Z_new,invR_Z_eta,data,Cov_post_new,S_R1,V_R1,sigmasqf); % log posterior of new proposal
    Lik(i)=logpznew;
    logpz=logconditionalZ(Z,emulator_Z,invR_Z_eta,data,Cov_post,S_R1,V_R1,sigmasqf);% log posterior of current state.
    %[pold2new,pnew2old]
    if Cov_post_new<restrict_emu
        if  logpnew2old==-inf % this means the probability of Z_new goes to Z is so small that it hardly come back, under this situation, new update is accepted immediatly.
            accepted_num=accepted_num+1;
            Z=Z_new;
            Grad_Z=Grad_Z_new;
            emulator_Z=emulator_Z_new;
            %Obs(i,:)=emulator_Z_new;
            Cov_post=Cov_post_new;
        else
            Pr=exp((logpznew-logpz+logpnew2old-logpold2new));
            %Pr=exp((logpznew-logpz));
            u=rand(1);
            %Obs(i,:)=emulator_Z;
            if Pr>=u %accept new proposal
                Z=Z_new;
                accepted_num=accepted_num+1;
                Grad_Z=Grad_Z_new;
                emulator_Z=emulator_Z_new;
                %Obs(i,:)=emulator_Z_new;
                Cov_post=Cov_post_new;
            end
        end
    end
    sigmasqf=conditionalsigmasq_f(invR_Z_eta,Z_new,data);
    chain(i,:)=[Z,eta,sigmasqf]; % new sample
    Err1(i)=(mean(abs(Z*data.KL'-True.f_input)));
    Err2(i)=(mean((Z*data.KL'-True.f_input).^2))^(1/2);
    Errinfty(i)=(max(abs(Z*data.KL'-True.f_input)));
    time_Z=toc;
    time_used=time_used+time_Z;
    if rem(i,500)==0
    str=[num2str(floor(n_simu/1000-i/1000)) 'K steps, ' num2str(floor(time_used*(n_simu-i)/i)), 's.',' P_a: ' num2str(floor(100*accepted_num/i))];
    waitbar(i/n_simu,h,str);
    end
end
close(h);
% reesult records every thing which may be used in the future.
result.likelihood=Lik;
result.data=data;
result.True=True;
result.para_emulator=para_emulator;
result.Z_MAP=MAP_est(chain((n_simu/6+1):50:end,1:data.M));% MAP estimation of KL coefficients
%result.obs=Obs;
result.chain=chain;%MCMC result
result.accepted_rate=accepted_num/n_simu; % acceptance rate of the Metropolis-Hastings algorithm for sampling from the full conditional distribution of \zeta
jump=floor((para_MCMC.nsimu-para_MCMC.burnin+1)/2500);
result.Z_postmean=mean(chain(para_MCMC.burnin:jump:end,1:data.M));% posterior mean of KL coefficients
result.f_postmean=data.KL*result.Z_postmean';% posterior mean of functional input
result.time=sum(time_used);% cpu time used for the whole MCMC
result.Var_pred=Var_pred; % recording of prediction variance of MCMC sample
result.Err1=Err1;%l_1 error of posterior sample
result.Err2=Err2;%l_2 error of posterior sample
result.Errinf=Errinfty;%l_infinity error of posterior sample
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
    u=rand(1)*conditionaleta(eta_old,Z,invR,logdetR,data,sigmasqf);
    L=eta_L;
    R=eta_U;
    temp=rand(1)*(R-L)+L;
    u_temp=conditionaleta(temp,Z,invR,logdetR,data,sigmasqf);
    while u_temp<u
        if temp<eta_old
           L=max(L,temp); 
        end
        if temp>eta_old
           R=min(R,temp);            
        end
        temp=rand(1)*(R-L)+L;
        u_temp=conditionaleta(temp,Z,invR,logdetR,data,sigmasqf);
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
function y=conditionaleta(eta,Z,invR,logdetR,data,sigmasqf)% using slise sampling
ind=round((eta-data.eta_L)*9800/(data.eta_U-data.eta_L)+1);
invSigma_z=invR{ind}/sigmasqf;
y=exp(-Z*(invSigma_z*Z')/2)*exp(-logdetR(ind)/2);
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
function y=conditionaleta_multi(eta,Z,invR,logdetR,data,sigmasqf)% using slise sampling
ind1=round((eta(1)-data.eta_L)*100+1);
ind2=round((eta(2)-data.eta_L)*100+1);
invSigma_z=invR{ind1,ind2}/sigmasqf;

y=exp(-Z*(invSigma_z*Z')/2)*(exp(logdetR(ind1,ind2)))^(-1/2); %likelihood
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
function y=logconditionalZ(Z,emulator_Z,invR_Z_eta,data,Cov_post,S_R1,V_R1,sigmasqf)
M=length(Z);
n=length(data.obs);
Sigma_e=Cov_post*diag(V_R1)+data.sigma_e^2;
logdetSig=sum(log(Sigma_e));
logprior=-(data.para_prior.a+M/2)*sigmasqf-(Z*(invR_Z_eta*Z')/(2*sigmasqf)+data.para_prior.b);
loglikelihood=-sum((data.obsp-emulator_Z).^2./Sigma_e')/2-logdetSig/2;
y=loglikelihood+logprior; %logposterior
end
function [y,sigma]=proposal(Z,invR_Z_eta,Grad,data,Cov_post,S_R1,V_R1,sigmasqf) % generating new proposal
global scale
n=length(data.obs);
Sigma_e=Cov_post*V_R1+data.sigma_e^2*eye(n);
sigma=invandlogdet(Grad*(Sigma_e\(Grad'))+invR_Z_eta/sigmasqf);
sigma=sigma/scale;
%mu=sigma*(Grad*(Sigma_e\(data.obs'-emulator_Z'+Grad'*Z')));
% global Mu
% Mu=[Mu;mu'];
y=mvnrnd(zeros(1,data.M),sigma)+Z;
end
function y=logtransform(Z_new,Z,invR_Z_eta,Grad,data,Cov_post,S_R1,V_R1,sigmasqf) 
global scale
n=length(data.obs);
Sigma_e=Cov_post*V_R1+data.sigma_e^2*eye(n);
[sigma,logdetinvsig]=invandlogdet(Grad*(Sigma_e\(Grad'))+invR_Z_eta/sigmasqf);
sigma=sigma/scale;
%mu=sigma*(Grad*(Sigma_e\(data.obs'-emulator_Z'+Grad'*Z')));
y=-(Z_new-Z)*((sigma)\(Z_new'-Z'))/2+logdetinvsig/2;
end
function y_grad=emulator_grad(x,para_emulator)
D_in=para_emulator.D_in;
D_para=para_emulator.D_para;
cor_func=para_emulator.cor_func;
thetaopt=para_emulator.thetaopt;
%R1=para_emulator.R1;
ResR2Ip=para_emulator.ResR2Ip;
%MTRIMI=para_emulator.MTRIMI;
[n1, d1]=size(D_in);
[n2, d2]=size(D_para);
R3=zeros(length(x),n2);
for j=1:n2
    R3(:,j)=CompCorr_diff(x,D_para(j,:),thetaopt((d1+1):end),cor_func);%(CompCorr(x_add,D_para(j,:),thetaopt((d1+1):end),cor_func)-CompCorr(x,D_para(j,:),thetaopt((d1+1):end),cor_func))/e;
end
y_grad=(ResR2Ip*(R3'))';
end
%%



function [y_pred,PosVar]=emulator_trans(x,para_emulator)
D_in=para_emulator.D_in;
D_para=para_emulator.D_para;
cor_func=para_emulator.cor_func;
thetaopt=para_emulator.thetaopt;
alphap=para_emulator.alphap;
sigma2=para_emulator.sigma2;
%R1=para_emulator.R1;
ResR2Ip=para_emulator.ResR2Ip;
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
Ypred=alphap+ResR2Ip*(R3');
y_pred=Ypred';
PosVar=max(sigma2*(1-R3*R2I*R3')+sigma2*(1-M2TRI*R3')^2/(sum(M2TRI)),0);
end
function y=CompCorr(x1,x2,eta,cor_func)
if(cor_func==0)
    y=prod(eta.^((x1-x2).^2));
elseif(cor_func==1)
    rho=x1-x2;
    rho1=sqrt(6)*abs(rho)./eta;
    y=prod((exp(-rho1)).*(rho1+1));
elseif(cor_func==2)
    rho=x1-x2;
    rho1=2*sqrt(2.5)*abs(rho)./eta;
    y=prod((exp(-rho1)/3).*(rho1.^2+3*rho1+3)); 
end
end
function y=CompCorr_1(x1,D,eta)
n=size(D,1);
    rho1=sqrt(6)*abs(x1-D)./(ones(n,1)*eta);
    y=prod((exp(-rho1')).*(rho1'+1));
end
function y=CompCorr_diff(x1,x2,eta,cor_func)
M=length(x1);
for i=1:M
if(cor_func==0)
    y=prod(eta.^((x1-x2).^2));
elseif(cor_func==1)
    rho=x1-x2;
    rho1=sqrt(6)*abs(rho)./eta;
    y(i)=(-6*(rho(i))/(eta(i)^2))*prod((exp(-rho1)).*(rho1+1))/(1+rho1(i));
elseif(cor_func==2)
    rho=x1-x2;
    rho1=2*sqrt(2.5)*abs(rho)./eta;
    y=prod((exp(-rho1)/3).*(rho1.^2+3*rho1+3)); 
end
end
y=y';
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