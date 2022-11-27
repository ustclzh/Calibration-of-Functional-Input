function [para_design,para_emulator,c_trace,Z_mode]=fully_sequential_generation(data,para_emulator,para_design,z_initial,criteria)
% d_desire: size of follow-up design size design size
% criteria: indicate the design criterion used in simulation, 1: WPV; 2: VL; 3: AEI;
%%
z_post=z_initial;
para_emulator.sumM2TRI=sum(para_emulator.M2TRI);
d_desire=para_design.d_desire;
R1=para_emulator.R1;
[s,V,d]=svd(R1);
para_emulator.s=s;
para_emulator.V=V;
data.obsp=data.obs*s;
para_emulator.alphap=s'*para_emulator.alpha;
para_emulator.ResR2Ip=s'*para_emulator.ResR2I;
Z_mode=[];
c_trace=[];
global optimal_design 
optimal_design.para_emulator=para_emulator;
optimal_design.data=data;
optimal_design.ss_min=0;
s=0;
while s<d_desire
    s=s+1;
    for i=1:size(para_emulator.D_para,1)
        [mk,~]=emulator_trans(para_emulator.D_para(i,:),para_emulator);
        mm(i)=sum(((data.obsp-mk).^2));
    end
    optimal_design.ss_min=min(mm);
    [z_add,para_emulator,z_post]=sequantially_add_one_point(data,para_emulator,z_post,criteria);
    [~,v]=emulator_trans(z_post,para_emulator);
    c_trace(s)=v;
    Z_mode=[Z_mode;z_post];
end
para_design.D_para=[para_emulator.D_para];
end

function [z_add,para_emulator,z_post]=sequantially_add_one_point(data,para_emulator,z_post,criteria)
%%
global optimal_design
optimal_design.para_emulator=para_emulator;
optimal_design.data=data;
switch criteria
    case 1
        if data.fix_eta==0
            if data.isotropic==0
                options = optimoptions('particleswarm','Display','off','InitialSwarmMatrix',[z_post],'HybridFcn',@patternsearch);
                [z_add,fval,exitflag] = particleswarm(@PV_criterion_aniso,data.M,[-10*ones(1,data.M)],[10*ones(1,data.M)],options);
                if max(abs(z_add))>9.9
                    options_patternsearch = optimoptions('patternsearch','Display','off','MaxIterations',2000*data.M);
                    [z_add1,fval2,exitflag] = patternsearch(@PV_criterion_aniso,[z_post+randn(1,data.M)/100],[],[],[],[],[-10*ones(1,data.M)],[10*ones(1,data.M)],[],options_patternsearch);
                    if fval2<fval
                        z_add=z_add1;
                    end
                end
            else
                 options = optimoptions('particleswarm','Display','off','InitialSwarmMatrix',[z_post],'HybridFcn',@patternsearch);
                [z_add,fval,exitflag] = particleswarm(@PV_criterion,data.M,[-10*ones(1,data.M)],[10*ones(1,data.M)],options);
                if max(abs(z_add))>9.9
                    options_patternsearch = optimoptions('patternsearch','Display','off','MaxIterations',2000*data.M);
                    [z_add1,fval2,exitflag] = patternsearch(@PV_criterion,[z_post+randn(1,data.M)/100],[],[],[],[],[-10*ones(1,data.M)],[10*ones(1,data.M)],[],options_patternsearch);
                    if fval2<fval
                        z_add=z_add1;
                    end
                end
            end
        else
            options = optimoptions('particleswarm','Display','off','InitialSwarmMatrix',[z_post],'HybridFcn',@patternsearch);
            [z_add,fval,exitflag] = particleswarm(@PV_criterion_fix_eta,data.M,[-10*ones(1,data.M)],[10*ones(1,data.M)],options);
            if max(abs(z_add))>9.9
                options_patternsearch = optimoptions('patternsearch','Display','off','MaxIterations',2000*data.M);
                [z_add1,fval2,exitflag] = patternsearch(@PV_criterion_fix_eta,[z_post+randn(1,data.M)/100],[],[],[],[],[-10*ones(1,data.M)],[10*ones(1,data.M)],[],options_patternsearch);
                if fval2<fval
                    z_add=z_add1;
                end
            end
        end
        z_add=z_add(1:data.M);
        y_add=Simulator(z_add,data);
        para_emulator=emulator_update(z_add,para_emulator,y_add);
    case 2
        options = optimoptions('particleswarm','Display','off','InitialSwarmMatrix',[z_post],'HybridFcn',@patternsearch);
        [z_add,fval,exitflag] = particleswarm(@VL_criterion,data.M,[-10*ones(1,data.M)],[10*ones(1,data.M)],options);
        if max(abs(z_add))>9.9
            options_patternsearch = optimoptions('patternsearch','Display','off','MaxIterations',2000*data.M);
            [z_add1,fval2,exitflag] = patternsearch(@VL_criterion,z_add,[],[],[],[],[-10*ones(1,data.M)],[10*ones(1,data.M)],[],options_patternsearch);
            if fval2<fval
                z_add=z_add1;
            end
        end
        z_add=z_add(1:data.M);
        y_add=Simulator(z_add,data);
        para_emulator=emulator_update(z_add,para_emulator,y_add);
    case 3
        options = optimoptions('particleswarm','Display','off','InitialSwarmMatrix',[z_post],'HybridFcn',@patternsearch);
        [z_add,fval,exitflag] = particleswarm(@AEI_criterion,data.M,[-10*ones(1,data.M)],[10*ones(1,data.M)],options);
        if max(abs(z_add))>9.9
            options_patternsearch = optimoptions('patternsearch','Display','off','MaxIterations',2000*data.M);
            [z_add1,fval2,exitflag] = patternsearch(@AEI_criterion,[z_post+randn(1,data.M)/100],[],[],[],[],[-10*ones(1,data.M)],[10*ones(1,data.M)],[],options_patternsearch);
            if fval2<fval
                z_add=z_add1;
            end
        end
        z_add=z_add(1:data.M);
        y_add=Simulator(z_add,data);
        para_emulator=emulator_update(z_add,para_emulator,y_add);
end
end
%%
function y=PV_criterion(para)
global optimal_design 
para_emulator=optimal_design.para_emulator;
data=optimal_design.data;
M=data.M;
yp=data.obsp;
sigma_e=data.sigma_e;
V=diag(para_emulator.V)';
para=para(1:M);
[mk,v]=emulator_trans(para,para_emulator);
logdets=sum(log(2*pi*(v*V+sigma_e^2)));
% invs2=data.invR_Z_Eta{ind};
% logs2=data.logdetR_Z_Eta(ind);
minuslogprior=-log(prior_marginal_z(para,data));%logs2/2+(M/2+data.para_prior.a)*log(para*invs2*para'/2+data.para_prior.b);
minusloglik=logdets/2+0.5*sum((yp-mk).^2./(v*V+sigma_e^2));
%y=-log(v)-minusloglik-minuslogprior;
y=-v*(exp(-minusloglik-minuslogprior));
end

function y=PV_criterion_aniso(para)
global optimal_design 
para_emulator=optimal_design.para_emulator;
data=optimal_design.data;
M=data.M;
yp=data.obsp;
sigma_e=data.sigma_e;
V=diag(para_emulator.V)';
para=para(1:M);
[mk,v]=emulator_trans(para,para_emulator);
logdets=sum(log(2*pi*(v*V+sigma_e^2)));
% invs2=data.invR_Z_Eta{ind};
% logs2=data.logdetR_Z_Eta(ind);
minuslogprior=-log(prior_marginal_z_aniso(para,data));%logs2/2+(M/2+data.para_prior.a)*log(para*invs2*para'/2+data.para_prior.b);
minusloglik=logdets/2+0.5*sum((yp-mk).^2./(v*V+sigma_e^2));
%y=-log(v)-minusloglik-minuslogprior;
y=-v*(exp(-minusloglik-minuslogprior));
end
function y=PV_criterion_fix_eta(para)
global optimal_design 
para_emulator=optimal_design.para_emulator;
data=optimal_design.data;
M=data.M;
yp=data.obsp;
sigma_e=data.sigma_e;
V=diag(para_emulator.V)';
para=para(1:M);
[mk,v]=emulator_trans(para,para_emulator);
logdets=sum(log(2*pi*(v*V+sigma_e^2)));
minuslogprior=(M/2+data.para_prior.a)*log(para*para'/2+data.para_prior.b);
minusloglik=logdets/2+0.5*sum((yp-mk).^2./(v*V+sigma_e^2));
%y=-log(v)-minusloglik-minuslogprior;
y=-v*(exp(-minusloglik-minuslogprior));
end
function prior=prior_marginal_z(z,data)
temp=z*(reshape(data.Inv_Sigma_eta'*z',data.M,99));
prior=mean(exp(-data.Log_Det_Sigma_eta/2).*((temp/2++data.para_prior.b).^(-data.M/2-data.para_prior.a)));
end
function prior=prior_marginal_z_aniso(z,data)
temp=z*(reshape(data.Inv_Sigma_eta'*z',data.M,201));
prior=mean(exp(-data.Log_Det_Sigma_eta/2).*((temp/2++data.para_prior.b).^(-data.M/2-data.para_prior.a)));
end
%%
function y=AEI_criterion(para)
global optimal_design
para_emulator=optimal_design.para_emulator;
data=optimal_design.data;
ss_min=optimal_design.ss_min;
yp=data.obsp;
n=length(yp);
R1=para_emulator.R1;
V=(diag(para_emulator.V))';
[mk,v]=emulator_trans(para,para_emulator);
sigmasq=2*v^2*sum(V.^2)+4*v*sum(((yp-mk).^2).*V);
mu=ss_min-n*v-sum((yp-mk).^2);
y=(-mu*normcdf(mu/sigmasq^(1/2))-sigmasq^(1/2)*normpdf(mu/sigmasq^(1/2)))*(v<1);
% y=-sigmasq^(1/2)*exp(-mu^2/sigmasq)/(2*pi)^(1/2)-mu*normcdf(mu/sigmasq^(1/2));
end
%%
function y=VL_criterion(para)
global optimal_design
para_emulator=optimal_design.para_emulator;
data=optimal_design.data;
yp=data.obsp;
n=length(yp);
sigma_e=data.sigma_e;
V=(diag(para_emulator.V))';
[mk,v]=emulator_trans(para,para_emulator);
dets1=prod(2*pi*(2*v*V+sigma_e^2));
dets2=prod(2*pi*(v*V+sigma_e^2));
y=((2*pi*sigma_e^2)^(-n/2)*(dets1)^(-1/2)*exp(-sum((yp-mk).^2./(2*v*V+sigma_e^2)))-(dets2)^(-1)*exp(-sum((yp-mk).^2./(v*V+sigma_e^2))));
y=-log(y+1);
end

%%
function [y_pred,PosVar,m2]=emulator_trans(x,para_emulator,anything)
D_in=para_emulator.D_in;
D_para=para_emulator.D_para;
%cor_func=para_emulator.cor_func;
thetaopt=para_emulator.thetaopt;
alphap=para_emulator.alphap;
sigma2=para_emulator.sigma2;
%R1=para_emulator.R1;
ResR2Ip=para_emulator.ResR2Ip;
R2I=para_emulator.R2I;
%MTRIMI=para_emulator.MTRIMI;
M2TRI=para_emulator.M2TRI;
sumM2TRI=para_emulator.sumM2TRI;
[n1, d1]=size(D_in);
%[n2, d2]=size(D_para);
% R3=zeros(1,n2);
% for j=1:n2
%     R3(j)=CompCorr(x,D_para(j,:),thetaopt((d1+1):end),cor_func);
% end
R3=CompCorr_1(x,D_para,thetaopt((d1+1):end));
y_pred=alphap'+R3*ResR2Ip';
PosVar=max(sigma2*(1-R3*R2I*R3')+sigma2*(1-M2TRI*R3')^2/(sumM2TRI),0);
if nargin==3
   m=para_emulator.m1+para_emulator.m2 *R3';
   m2=m(end);
end
end

function [y_pred,PosVar,m2]=emulator(x,para_emulator,anything)
D_in=para_emulator.D_in;
D_para=para_emulator.D_para;
%cor_func=para_emulator.cor_func;
thetaopt=para_emulator.thetaopt;
alpha=para_emulator.alpha;
sigma2=para_emulator.sigma2;
%R1=para_emulator.R1;
ResR2I=para_emulator.ResR2I;
R2I=para_emulator.R2I;
%MTRIMI=para_emulator.MTRIMI;
M2TRI=para_emulator.M2TRI;
sumM2TRI=para_emulator.sumM2TRI;
[n1, d1]=size(D_in);
%[n2, d2]=size(D_para);
% R3=zeros(1,n2);
% for j=1:n2
%     R3(j)=CompCorr(x,D_para(j,:),thetaopt((d1+1):end),cor_func);
% end
R3=CompCorr_1(x,D_para,thetaopt((d1+1):end));
y_pred=alpha'+R3*ResR2I';
PosVar=max(sigma2*(1-R3*R2I*R3')+sigma2*(1-M2TRI*R3')^2/(sumM2TRI),0);
if nargin==3
   m=para_emulator_add.m1+para_emulator.m2 *R3';
   m2=m(end);
end
end
%%
function para_emulator=emulator_update(z_add,para_emulator,y_add)
D_para=para_emulator.D_para;
d1=size(para_emulator.D_in,2);
d2=size(para_emulator.D_para,2);
%R3=CompCorr_1(z_add,D_para,para_emulator.thetaopt((d1+1):(d1+d2)));
D_para=[D_para;z_add];
Y=para_emulator.Y;
Y=[Y,y_add'];
n2=size(D_para,1)-1;
R3=zeros(1,n2);
n1=size(para_emulator.D_in,1);
for j=1:n2
    R3(j)=CompCorr(z_add,D_para(j,:),para_emulator.thetaopt((d1+1):end),para_emulator.cor_func);
end
%R2I=block_matrix_inversion(para_emulator.R2I,R3,1);
R2=[para_emulator.R2,R3';[R3,1]];
R2I=invandlogdet(R2);
R1=para_emulator.R1;
R1I=invandlogdet(R1);
n2=size(D_para,1);
M2TRI=ones(1,n2)*R2I;
SSR2I=M2TRI*ones(n2,1);
alpha=Y*(M2TRI')/SSR2I;
MTRIMI=R1/SSR2I;
Res=Y-alpha*ones(1,n2);
para_emulator.m1=(M2TRI')/SSR2I;
para_emulator.m2=(eye(n2)-((M2TRI')/SSR2I)*ones(1,n2))*R2I;

A=Res'*R1I;
ResR2I=Res*R2I;
sigma2=sum(sum(A'.*ResR2I))/(n1*n2-n2);
para_emulator.D_para=D_para;
para_emulator.alpha=alpha;
para_emulator.sigma2=sigma2;
para_emulator.R1=R1;
para_emulator.R2=R2;
para_emulator.ResR2I=ResR2I;
para_emulator.R2I=R2I;
para_emulator.MTRIMI=MTRIMI;
para_emulator.M2TRI=M2TRI;
para_emulator.sumM2TRI=sum(M2TRI);
para_emulator.Y=Y;
[s,V,d]=svd(R1);
para_emulator.alphap=s'*para_emulator.alpha;
para_emulator.ResR2Ip=s'*para_emulator.ResR2I;
end
%%
function y=CompCorr(x1,x2,eta,cor_func)
global nugterm
if(cor_func==0)
    y=prod(eta.^((x1-x2).^2));
elseif(cor_func==1)
    rho=x1-x2;
    rho1=sqrt(6)*abs(rho)./eta;
    y=prod((exp(-rho1)).*(rho1+1))+nugterm*prod(x1==x2);
elseif(cor_func==2)
    rho=x1-x2;
    rho1=2*sqrt(2.5)*abs(rho)./eta;
    y=prod((exp(-rho1)/3).*(rho1.^2+3*rho1+3)); 
end
end

function y=CompCorr_1(x1,D,eta)
global nugterm
    rho1=sqrt(6)*abs(x1-D)./eta;
    y=prod((exp(-rho1')).*(rho1'+1))+nugterm*prod((x1==D)');
end

