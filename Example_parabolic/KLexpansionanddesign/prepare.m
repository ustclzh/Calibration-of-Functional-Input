
M=9;
sigma_e=0.02;
iso=1;
%pre_KL=preparition_for_general_data(sigma_e,M,1,iso);
pre_KLL=preparition_for_general_data(sigma_e,M,2,iso);
pre_KLC=preparition_for_general_data(sigma_e,M,3,iso);
pre_KLU=preparition_for_general_data(sigma_e,M,4,iso);
iso=0;
pre_KL_aniso=preparition_for_general_data(sigma_e,M,5,iso);
%filename=['pre', num2str(M),'basis'];
%save(filename,'-v7.3')
data1.M=50;
d=randn(3000,data1.M);
para_design.d_desire=300;
para_design.cor_func=1;
para_design.initial_design=[];
para_design.candidate_sample=d;
data1.D_in=[];
para_design=Generate_Design(data1,para_design,'sp');
design_z=para_design.D_para;
%save('inputdesign.mat','design_z','-v7.3')

function y=preparition_for_general_data(sigma_e,M,basis,iso)
[True,data]=Pde_parameter(0.5,zeros(1,100),sigma_e,iso);
data=generate_basis_all(data,basis,M);
multi_design=10;
para_MCMC.nsimu=60000;
para_MCMC.fixeta=data.fix_eta;
if data.fix_eta==1
    para_MCMC.eta=data.eta;
end
candidate_sample=Sample_from_prior(data,para_MCMC);%2500 candidate
para_design.d_desire=data.M*multi_design;
para_design.cor_func=1;
para_design.initial_design=[];
para_design.candidate_sample=candidate_sample;
para_design=Generate_Design(data,para_design,'fr');
para_emulator=Train_Emulator(data,para_design);
data.obsp=data.obs*para_emulator.s;
y.data=data;
y.para_emulator=para_emulator;
y.Candidate_for_start=candidate_sample;
end
%%
function chain=Sample_from_prior(data,para_MCMC)
% para_MCMC.nsimu: length of MCMC
% start_point: starting point of the chain
% para_MCMC.fixeta: an indicator for whether fixing eta. 1 for fixed 0 for
% not fixed
% para_MCMC.eta : the eta value that is fixed. Only valid when fixeta=1.
%% preparation
time_used=0;
n_simu=para_MCMC.nsimu;
M=data.M;
p=M+3-data.isotropic;
if data.fix_eta==1
   eta=data.eta;
   R_Z_eta=eye(M);
   [S,V,D]=svd(R_Z_eta);
   R_z_trans=S*V^(1/2)*S';
else 
    eta=0.5*ones(1,2-data.isotropic)*(data.eta_U-data.eta_L)+data.eta_L;
    invR_Z_Eta=data.invR_Z_Eta;
    logdetR_Z_Eta=data.logdetR_Z_Eta;
    R_Z_Eta=data.R_Z_Eta;
end
chain=zeros(n_simu,p);
sigmasq_f=1;
Z=randn(1,M);
for i=1:n_simu
    tic
    if data.fix_eta==0
        if data.isotropic==1
            eta=slicesampling(eta,sigmasq_f,Z,invR_Z_Eta,logdetR_Z_Eta,data);
            ind=floor((eta-data.eta_L)*9800/(data.eta_U-data.eta_L)+1.0001);
            R_Z_eta=R_Z_Eta{ind};
            [S,V,D]=svd(R_Z_eta);
            R_z_trans=S*V^(1/2)*S';
            sigmasq_f=conditionalsigmasq_f(R_Z_eta,Z,data); % update sigmasq_f
        else
            eta=slicesampling_multi(eta,sigmasq_f,Z,invR_Z_Eta,logdetR_Z_Eta,data);
            ind1=round((eta(1)-data.eta_L)*98/(data.eta_U-data.eta_L)+1);
            ind2=round((eta(2)-data.eta_L)*98/(data.eta_U-data.eta_L)+1);
            R_Z_eta=R_Z_Eta{ind1,ind2};
            [S,V,D]=svd(R_Z_eta);
            R_z_trans=S*V^(1/2)*S';
            sigmasq_f=conditionalsigmasq_f(R_Z_eta,Z,data); % update sigmasq_f
        end
    else
        sigmasq_f=conditionalsigmasq_f(R_Z_eta,Z,data); % update sigmasq_f
    end
    temp=randn(1,data.M);
    Z_new=(sigmasq_f^(1/2)*R_z_trans*temp')';
    chain(i,:)=[Z_new,sigmasq_f,eta]; % new sample
    time_Z=toc;
    time_used=time_used+time_Z;
end
chain=chain(10001:20:end,1:M);
end
%slice sampling
function y=slicesampling(eta_old,sigmasq_f,Z,invR,logdetR,data) 
eta_L=data.eta_L;
eta_U=data.eta_U;
% slice sampling 
   for i=1:10
    u=rand(1)*conditionaleta(eta_old,sigmasq_f,Z,invR,logdetR,data);
    L=eta_L;
    R=eta_U;
    temp=rand(1)*(R-L)+L;
    u_temp=conditionaleta(temp,sigmasq_f,Z,invR,logdetR,data);
    while u_temp<u
        if temp<eta_old
           L=max(L,temp); 
        end
        if temp>eta_old
           R=min(R,temp);            
        end
        temp=rand(1)*(R-L)+L;
        u_temp=conditionaleta(temp,sigmasq_f,Z,invR,logdetR,data);
    end
    eta_new=temp;
    eta_old=eta_new;
   end
    if eta_new<eta_L
        eta_new=eta_L;
    elseif eta_new>eta_U
        eta_new=eta_U;
    end
    y=eta_new;
end
function y=conditionaleta(eta,sigmasq_f,Z,invR,logdetR,data)
ind=round((eta-data.eta_L)*9800/(data.eta_U-data.eta_L)+1);
invSigma_z=invR{ind}/sigmasq_f;
y=exp(-Z*(invSigma_z*Z')/2)*(exp(logdetR(ind)))^(-1/2); %likelihood
end
%
function y=slicesampling_multi(eta_old,sigmasq_f,Z,invR,logdetR,data) 
eta_L=data.eta_L;
eta_U=data.eta_U;
% slice sampling 
for j=1:2
    eta_temp=eta_old;
   for i=1:10
    u=rand(1)*conditionaleta_multi(eta_old,sigmasq_f,Z,invR,logdetR,data);
    L=eta_L;
    R=eta_U;
    temp=rand(1)*(R-L)+L;
    eta_temp(j)=temp;
    u_temp=conditionaleta_multi(eta_temp,sigmasq_f,Z,invR,logdetR,data);
    while u_temp<u
        if temp<eta_old(j)
           L=max(L,temp); 
        end
        if temp>eta_old(j)
           R=min(R,temp);            
        end
        temp=rand(1)*(R-L)+L;
        eta_temp(j)=temp;
        u_temp=conditionaleta_multi(eta_temp,sigmasq_f,Z,invR,logdetR,data);
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
function y=conditionaleta_multi(eta,sigmasq_f,Z,invR,logdetR,data)% using slise sampling
ind1=round((eta(1)-data.eta_L)*98/(data.eta_U-data.eta_L)+1);
ind2=round((eta(2)-data.eta_L)*98/(data.eta_U-data.eta_L)+1);
invSigma_z=invR{ind1,ind2}/sigmasq_f;
y=exp(-Z*(invSigma_z*Z')/2)*(exp(logdetR(ind1,ind2)))^(-1/2); %likelihood
end
%
function y=conditionalsigmasq_f(R_zeta_eta,Z,data) 
%return an inverse gamma sample
% the prior is invgamma(a,b)
a=data.para_prior.a; 
b=data.para_prior.b; 
a=a+data.M/2;
b=b+Z*(R_zeta_eta\Z')/2;
y=(gamrnd(a,1/b))^(-1);
if y>15
    y=15;
end
if isnan(y)
    y=15;
end
end
%%
function para_design=Generate_Design(data,para_design,algo)
% check
% data:necessary data;
% para_design.candidate_sample: candidate sample for generating design
% para_emulator: all parameters for emulator
% follow up design;
% design 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if algo=='fr'% minimax design for finite region
% if algo=='sp'% support points 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
candidate_sample=para_design.candidate_sample(:,1:data.M);
if size(candidate_sample,1)>10000
     candidate_sample=candidate_sample(end-10000:end,:);
end
candidate_sample=dupliction_removing(candidate_sample); % remove the same sampel in candidate.
D_in=data.D_in;
d_desire=para_design.d_desire;
d_desire=max(d_desire,50);
MaxIter=200;
if algo=='fr'% minimax design for finite region
    D_para=generate_design_finite_region(candidate_sample,d_desire);
    if size(D_para,1)>d_desire
        D_para=D_para(end-d_desire+1:end,:);
    end
    D_para=[para_design.initial_design;D_para];
end
if algo=='sp'%support points 
    D_para=generate_design_support_points(d_desire,candidate_sample,MaxIter);
    D_para=[para_design.initial_design;D_para];
end
if isempty(D_para)
    disp('Generate design failed, using initial sample directly');
    D_para=candidate_sample;
    if size(D_para,1)>20*data.M
        D_para=D_para(1:20*data.M,:);
    end
end
para_design.D_in=D_in;
para_design.D_para=D_para;
end
%
function D_candidate=dupliction_removing(sample)
N=size(sample,1);
Dist_matirx=zeros(N,N);
for i=1:N
    for j=1:N
        Dist_matirx(i,j)=sum((sample(i,:)-sample(j,:)).^2);
    end
end
Index=zeros(1,N);
for i=1:N-1
    for j=i+1:N
        if Dist_matirx(i,j)<0.0001
            Index(j)=Index(j)+1;
        end
    end
end
D_candidate=sample(Index==0,:);
end
%support point
function y=generate_design_support_points(n,Sample,maxIte)
% check
% warning: chack that Sample doesn't contain any repeatation
% support points
% n: desired sample size
% Sample: candidate sample, i.e. MC sample of a target distribution
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
%finite region
function D_para=generate_design_finite_region(sample,d_desire)
% check 1
% design generated from prior sample
% sample: candiate sample, sample size 2000;
% d_desire: desired design size
if size(sample,1)>2000
sample=sample(end-2000+1:end,:);
end
%%
Design_index=Minimax_Design_from_Sample(sample,d_desire,0.6,1.2,2,2);
D_para=sample(Design_index,:);
end
function NearMinimaxSolution=Minimax_Design_from_Sample(CandidateSet,n,fl,fu,reduce,DP)
%ProcedureC(CandidateSet,n,fl,fu,reduce,DP)
%This program implements a heuristic procedure that finds near-minimax designs.
%CandidateSet should be a matlab variable that contains the candidate set.
%n=desired design size
%L=fl*DKS, U=fu*DKS, DKS=distance of n-point KS design 
%If reduce=1, then reduced SCLP is solved. Otherwise, SCLP is solved.
%All distances are computed to DP decimal places.
global N_fr StandardizedCS Distmat_fr 
format long g
if(n<=2||(ismember(mod(n,2),[0 1])~=1))
    disp('error')
    return
end
[N_fr, ~]=size(CandidateSet);
StandardizedCS=(CandidateSet-repmat(mean(CandidateSet),N_fr,1))./repmat(std(CandidateSet,1,1),N_fr,1);
Distmat_fr=pdist2(StandardizedCS,StandardizedCS);
Distmat_fr=round(Distmat_fr*10^DP)/10^DP;
[MaxDistforEachCol, RowIndices]=max(Distmat_fr);
[~, MaxDistColumnIndex]=max(MaxDistforEachCol);
MaxDistPair=[MaxDistColumnIndex RowIndices(MaxDistColumnIndex)];
[~, InitialDesignIndices]=SequentialHeuristicDesign2(MaxDistPair,n);
InitialNonDesignIndices=setdiff(1:N_fr,InitialDesignIndices);
InitialDist=max(min(Distmat_fr(InitialDesignIndices,InitialNonDesignIndices)));
U=fu*InitialDist;
L=fl*InitialDist ;
Phip=unique(Distmat_fr);
m=length(Phip);
theta=zeros(m,1);
for i=1:(m-1)
    theta(i)=(Phip(i)+Phip(i+1))/2;
end    
    theta(m)=1.1*Phip(m);
index=(L<=Phip)&(Phip<=U);
Sval=sort(theta(index),'descend');
LS=length(Sval);
Candidatedist=sort(Phip(index),'descend');
storexopt=zeros(LS,N_fr);
storez=zeros(LS,1);
for i=1:LS
if(i==1)
    check=1;
else
    check=designdist>Sval(i);
end
if(check)
    indexmat=Distmat_fr<=Sval(i);
    if(reduce==1)
    [A, remainingvar]=trimmatrix(indexmat,N_fr);
    A=logical(A);
    else
    A=logical(indexmat);
    end
    xopt=HGW(A);
    z=sum(xopt);
end
storez(i)=z;
    if(reduce==1)
    storexopt(i,remainingvar)=xopt;  
    else
    storexopt(i,:)=xopt;
    end
if(check)
designdist=max(min(Distmat_fr(logical(storexopt(i,:)),logical(storexopt(i,:)==0))));
end    
end
for i=1:LS-1
    for j=i+1:LS
        if(storez(i)>storez(j))
            storez(i)=storez(j);
        end
    end 
end
MaxDesignSizeIndex=storez==max(storez);
storez(MaxDesignSizeIndex)=[];
Candidatedist(MaxDesignSizeIndex)=[];
MinimaxDesignSizes=unique(storez);
NoDesignSizes=length(MinimaxDesignSizes);
MinimaxDistances=zeros(NoDesignSizes,1);
for i=1:NoDesignSizes
    IndexMinimaxDistance=find(abs(storez-MinimaxDesignSizes(i))<10^-5, 1, 'last' );
    MinimaxDistances(i)=Candidatedist(IndexMinimaxDistance);
end
DesignSize=find(MinimaxDesignSizes>=n,1);
DesignSize=MinimaxDesignSizes(DesignSize);
if isempty(DesignSize)
    DesignSize=MinimaxDesignSizes(end);
end
IndexNearMinimaxSolution=find(storez==DesignSize, 1, 'last' );    
NearMinimaxSolution=storexopt(IndexNearMinimaxSolution,:)>(1-10^-5);
end
function [Design, DesignIndices]=SequentialHeuristicDesign2(StartDesignIndices,n)
global StandardizedCS
n0=length(StartDesignIndices);
DesignIndices=StartDesignIndices;
subset=StandardizedCS(StartDesignIndices,:);
for i=1:(n-n0)
    NewPointIndex=findpoint(DesignIndices);
    DesignIndices=[DesignIndices NewPointIndex];
    subset=[subset; StandardizedCS(NewPointIndex,:)];
end 
Design=subset;
end   
function findpoint=findpoint(DesignIndices)
global Distmat_fr N_fr
NonDesignIndices=setdiff(1:N_fr,DesignIndices);
Dist=min(Distmat_fr(DesignIndices,NonDesignIndices));
[~, index]=max(Dist);
findpoint=NonDesignIndices(index);
end
function HGW=HGW(A)
[m, N]=size(A);
b=-ones(m,1);
c=ones(1,N);
options=optimset('LargeScale','on','Display','off');
[xopt , ~ , ~] = linprog(c,-A,b,[],[],zeros(1,N),ones(1,N),[],options);
v=max(sum(A,2));
tentative=xopt>=(1/v);
Xi=find(tentative);
n1=sum(tentative);
stop=0;
while(stop==0)
Anew=A(:,Xi);
redun=sum(Anew,2);
mr=zeros(n1,1);
for i=1:n1
    mr(i)=min(redun(Anew(:,i)));
end
[value, jstar]=max(mr);
if(value>1)
    tentative(Xi(jstar))=0;
    Xi(jstar)=[];
    n1=n1-1;
else
    stop=1;
end
end
HGW=tentative;
end
function [trimmedmatrix, remainingvar]=trimmatrix(indexmat,N)
i=1;
n=N;
while(i<n)
    check=0;
    j=i+1;
    while((check==0)&&(j<=n))
        Diff=indexmat(i,:)-indexmat(j,:);
        if(all(Diff>=0)&&(sum(Diff)>10^-6))
            indexmat(i,:)=[];
            n=n-1;
            check=1;
        elseif(all(Diff<=0))
            indexmat(j,:)=[];
            n=n-1;      
        else
            j=j+1;
        end
    end
    if(check==0)
        i=i+1;
    end
end

storeindex=1:N;
i=1;
n=N;
while(i<n)
    check=0;
    j=i+1;
    while((check==0)&&(j<=n))
        Diff=indexmat(:,j)-indexmat(:,i);
        if(all(Diff>=0)&&(sum(Diff)>10^-6))
            indexmat(:,i)=[];
            storeindex(i)=[];
            n=n-1;
            check=1;
        elseif(all(Diff<=0))
            indexmat(:,j)=[];
            storeindex(j)=[];
            n=n-1;                
        else
            j=j+1;
        end
    end
    if(check==0)
        i=i+1;
    end
end

trimmedmatrix=indexmat;
remainingvar=storeindex;
end




