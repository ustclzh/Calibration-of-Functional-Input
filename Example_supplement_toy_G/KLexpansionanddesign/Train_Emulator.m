function [para_emulator,data]=Train_Emulator(data,para_design,para_emulator)
% check
%% training the emulation
% para_design.D_in: design for space or space-time domain, i.e.,
% observation nodes
% para_design.D_para: design for KL coefficients
% para_design.cor_func: class of correlation funciton
%tic 
if nargin==3 % means update emulator
   Y_0=para_emulator.Y;
   q_0=size(Y_0,2);
   Y=simulator_out(para_design.D_para(q_0+1:end,:),data);
   Y=[Y_0,Y];
elseif nargin==2 % means first train
    %%
Y=simulator_out(para_design.D_para,data);
else
para_design=data.para_design;
Y=data.para_emulator.Y;
end
%time= toc;
% disp(['Run simulators on design spends ' num2str(time) 'seconds']);
%tic
%%
para_emulator=train_emulator(para_design.D_in,para_design.D_para,Y,para_design.cor_func);
%time= toc;
data.obsp=data.obs*para_emulator.s;
% disp(['Tran emulator time ' num2str(time) 'seconds']);
end
function Y=simulator_out(D_para,data)
%%
d=size(D_para,1);
n=length(data.obs);
Y=zeros(n,d);
for i=1:d
    Y(:,i)=Simulator(D_para(i,:),data)';
end
end
function para_emulator=train_emulator(D_in,D_para,Y,cor_func,start)
%%
global Y_simulator sample_all parameter_all Type_correlation_function sample_obs sample_parameter correlation_obs correlation_para scale_train Correlation_Matrix_obs Correlation_Matrix_para
Type_correlation_function=cor_func;
Y_simulator=Y;
[sample_obs, correlation_obs]=size(D_in);
[sample_parameter, correlation_para]=size(D_para);
sample_all=sample_obs*sample_parameter;
parameter_all=correlation_obs+correlation_para;
if(sum(abs(size(Y_simulator)-[sample_obs sample_parameter]))>10^-6)
    display('error')
    return    
end
scale_train=[range(D_in,1) range(D_para,1)];

Correlation_Matrix_obs=cell(1,correlation_obs);
for i=1:correlation_obs
    Correlation_Matrix_obs{i}=abs(repmat(D_in(:,i),1,sample_obs)-repmat(D_in(:,i)',sample_obs,1));
end
Correlation_Matrix_para=cell(1,correlation_para);
for i=1:correlation_para
    Correlation_Matrix_para{i}=abs(repmat(D_para(:,i),1,sample_parameter)-repmat(D_para(:,i)',sample_parameter,1));
end

options=optimoptions(@fminunc,'MaxIter',10^3,'TolX',0.01,'TolFun',0.01,'MaxFunEvals',10^5,'Display','off','Algorithm','quasi-newton','ObjectiveLimit',-10^250);

if(nargin==5)
    [paropt1, fval1, exitflag1]=fminunc(@Obj,start,options); 
    paropt2=paropt1;
    fval2=fval1;
    exitflag2=exitflag1;
else
if(Type_correlation_function==0)
    [paropt1, fval1, exitflag1]=fminunc(@Obj,[norminv(0.25.^(1./(scale_train.^2)))],options); 
    [paropt2, fval2, exitflag2]=fminunc(@Obj,[norminv(0.75.^(1./(scale_train.^2)))],options);
elseif(Type_correlation_function==1)
    [paropt1, fval1, exitflag1]=fminunc(@Obj,[log(0.909700041540068.*scale_train)],options);
    [paropt2, fval2, exitflag2]=fminunc(@Obj,[log(2.54815755509488.*scale_train)],options);  
elseif(Type_correlation_function==2)
    [paropt1, fval1, exitflag1]=fminunc(@Obj,log(0.888569182612929.*scale_train),options); 
    [paropt2, fval2, exitflag2]=fminunc(@Obj,log(2.24172067901287.*scale_train),options);    
end
end
if((exitflag1<=0)||(exitflag2<=0))
    display('Likelihood optimization failure: algorithm did not converge.')   
end
if(fval1<fval2)
    if(Type_correlation_function==0)
        thetaopt=normcdf(paropt1);
    elseif(Type_correlation_function==1)
        thetaopt=exp(paropt1);        
    elseif(Type_correlation_function==2)
        thetaopt=exp(paropt1);        
    end
else
    if(Type_correlation_function==0)
        thetaopt=normcdf(paropt2);
    elseif(Type_correlation_function==1)
        thetaopt=exp(paropt2);  
    elseif(Type_correlation_function==2)
        thetaopt=exp(paropt2);          
    end
end
% display('Optimal values of GP parameters')
% disp(thetaopt)
R1=CompCorr1(thetaopt(1:correlation_obs));
%R1=eye(size(R1,1));
R2=CompCorr2(thetaopt((correlation_obs+1):parameter_all));
R1I=invandlogdet(R1);
R2I=invandlogdet(R2);
M2TRI=ones(1,sample_parameter)*R2I;
SSR2I=M2TRI*ones(sample_parameter,1);
alpha=Y_simulator*(M2TRI')/SSR2I;
MTRIMI=R1/SSR2I;
Res=Y_simulator-alpha*ones(1,sample_parameter);
A=Res'*R1I;
ResR2I=Res*R2I;
sigma2=sum(sum(A'.*ResR2I))/(sample_all-sample_parameter);
[s,V,~]=svd(R1);
para_emulator.s=s;
para_emulator.V=V;
para_emulator.alphap=s'*alpha;
para_emulator.ResR2Ip=s'*ResR2I;
para_emulator.R2=R2;
para_emulator.D_in=D_in;
para_emulator.D_para=D_para;
para_emulator.cor_func=cor_func;
para_emulator.thetaopt=thetaopt;
para_emulator.alpha=alpha;
para_emulator.sigma2=sigma2;
para_emulator.R1=R1;
para_emulator.ResR2I=ResR2I;
para_emulator.R2I=R2I;
para_emulator.MTRIMI=MTRIMI;
para_emulator.M2TRI=M2TRI;
para_emulator.sumM2TRI=sum(M2TRI);
para_emulator.Y=Y;
end
function Objective=Obj(par)
global Y_simulator sample_all parameter_all Type_correlation_function sample_obs sample_parameter correlation_obs scale_train
if(sum(isnan(par))>0)
    Objective=Inf;
    return
end
if(Type_correlation_function==0)
    r=normcdf(par);
    if(sum(r>(0.999.^(1./(scale_train.^2))))>0)
        Objective=Inf;
        return
    end
elseif(Type_correlation_function==1)
    r=exp(par);    
    if(sum(r>(53.9511207457682.*scale_train))>0)
        Objective=Inf;
        return        
    end      
elseif(Type_correlation_function==2)
    r=exp(par);
    if(sum(r>(40.7953930912638.*scale_train))>0)
        Objective=Inf;
        return        
    end
end

R1=CompCorr1(r(1:correlation_obs));
R2=CompCorr2(r((correlation_obs+1):parameter_all));

[R1I,LD1]=invandlogdet(R1);
[R2I,LD2]=invandlogdet(R2);
LDR=(sample_parameter)*LD1+(sample_obs)*LD2;
SSR2I=sum(sum(R2I));
alpha=Y_simulator*R2I*ones(sample_parameter,1)/SSR2I;

Res=Y_simulator-alpha*ones(1,sample_parameter);
A=Res'*R1I;
B=Res*R2I;
sigma2=sum(sum(A'.*B))/(sample_all-sample_obs);
Objective=(sample_all-sample_obs)*log(max(sigma2,0))+LDR+sample_obs*log(SSR2I)-LD1;
if(Objective<(-10^250))
    Objective=-10^250;
elseif(isnan(Objective)||(Objective>10^250))
    Objective=10^250;    
end
end
function Corr=CompCorr1(r)
global Type_correlation_function Correlation_Matrix_obs correlation_obs

if(Type_correlation_function==0)
    Corr=1;
    for i=1:correlation_obs
        Corr=Corr.*(r(i).^((Correlation_Matrix_obs{i}).^2));
    end
elseif(Type_correlation_function==1)
	Corr=1;
    for i=1:correlation_obs
        rho1=sqrt(6)*Correlation_Matrix_obs{i}./r(i);
        Corr=Corr.*((exp(-rho1)).*(rho1+1));
    end
elseif(Type_correlation_function==2)
    Corr=1;
    for i=1:correlation_obs
        rho1=2*sqrt(2.5)*Correlation_Matrix_obs{i}./r(i);
        Corr=Corr.*((exp(-rho1)/3).*(rho1.^2+3*rho1+3)); 
    end
end
end
function Corr=CompCorr2(r)
global Type_correlation_function Correlation_Matrix_para correlation_para nugterm

if(Type_correlation_function==0)
    Corr=1;
    for i=1:correlation_para
        Corr=Corr.*(r(i).^((Correlation_Matrix_para{i}).^2));
    end
elseif(Type_correlation_function==1)
	Corr=1;
    for i=1:correlation_para
        rho1=sqrt(6)*Correlation_Matrix_para{i}./r(i);
        Corr=Corr.*((exp(-rho1)).*(rho1+1));
    end
    Corr=Corr+nugterm*eye(size(rho1,1));
elseif(Type_correlation_function==2)
    Corr=1;
    for i=1:correlation_para
        rho1=2*sqrt(2.5)*Correlation_Matrix_para{i}./r(i);
        Corr=Corr.*((exp(-rho1)/3).*(rho1.^2+3*rho1+3)); 
    end
end
end