function data=generate_basis_all(data,basis_type,M,eta,eta_prior)
% check
% return the struct data type
% parameter_setting: parameter including the PDE information, MCMC
% information, prior information
para_prior.a=3;
para_prior.b=3;
data.para_prior=para_prior;
N=length(data.X);
data.M=M;
switch basis_type
    case 1
        eta_all=[0:0.01:1]*(data.eta_U-data.eta_L)+data.eta_L;
        Sigma=zeros(N,N);
        R=data.R;
        for i=1:length(eta_all)
            Sigma=Sigma+R(data.X,data.Y,eta_all(i));
        end
        Sigma=Sigma/length(eta_all);
        [S,V,~]=svd(Sigma);
%         if M==size(data.X,1)
%             M=find(diag(V)<10^-5,1);
%             data.M=M;
%         end
%         M
        KL=S(:,1:M)*(V(1:M,1:M).^0.5);%each column is the principal components. 
        KL_trans=S(:,1:M)*(inv(V(1:M,1:M)).^0.5);
        %figure
        %plot(KL)
        KL_full=S*(V.^0.5);
        data.KL_full=KL_full;
        data.KL_trans=KL_trans;
        data.KL=KL;
        data.eig=diag(V);
        data.M=M;
        data.N=N;
        data.eigen=V;
        data.fix_eta=0;
    case 2 
        data=KLbasislower(data,M);
        data.fix_eta=1;
        data.eta=data.eta_L;
    case 3
        data=KLbasiscentral(data,M);
        data.fix_eta=1;
        data.eta=(data.eta_L+data.eta_U)/2;
    case 4
        data=KLbasisupper(data,M);  
        data.fix_eta=1;
        data.eta=data.eta_U;
    case 5 %anisotropic GP input
        temp=sobolset(2);
        temp=net(temp,1000);
        eta_all=temp*(data.eta_U-data.eta_L)+data.eta_L;
        Sigma=zeros(N,N);
        R=data.R;
        for i=1:1000
            Sigma=Sigma+R(data.X,data.Y,eta_all(i,:));
        end
        Sigma=Sigma/1000;
        [S,V,~]=svd(Sigma);
        KL=S(:,1:M)*(V(1:M,1:M).^0.5);%each column is the principal components. 
        KL_full=S*(V.^0.5);
        data.KL_full=KL_full;
        KL_trans=S(:,1:M)*(inv(V(1:M,1:M)).^0.5);
        data.KL_trans=KL_trans;
        data.KL=KL;
        data.eig=diag(V);
        data.M=M;
        data.N=N;
        data.eigen=V;
        data.fix_eta=0;
        a=sobolset(2);
        b=net(a,201);
        c=round(b*98)+1;
        data.ind=c;
        data.isotropic=0;

    case 6
        if nargin ==4
            data=KLbasiscustmize(data,M,eta);  
        else
            data=KLbasiscustmize(data,M,[0.99,0.01]);  
        end
        data.fix_eta=1;
        data.eta=data.eta_U;
    case 7 %other kind of prior
        if nargin==5
            eta_all=eta_prior;
        else
            eta_all=[0:0.01:1]*(data.eta_U-data.eta_L)+data.eta_L;
        end
        Sigma=zeros(N,N);
        R=data.R;
        for i=1:length(eta_all)
            Sigma=Sigma+R(data.X,data.Y,eta_all(i));
        end
        Sigma=Sigma/length(eta_all);
        [S,V,~]=svd(Sigma);
        KL=S(:,1:M)*(V(1:M,1:M).^0.5);%each column is the principal components. 
        KL_trans=S(:,1:M)*(inv(V(1:M,1:M)).^0.5);
        KL_full=S*(V.^0.5);
        data.KL_full=KL_full;
        data.KL_trans=KL_trans;
        data.KL=KL;
        data.eig=diag(V);
        data.M=M;
        data.N=N;
        data.eigen=V;
        data.fix_eta=0;
end  
if basis_type==1
    k=9801;
    invR_Z_Eta=cell(1,k);
    logdetR_Z_Eta=zeros(1,k);
    R_Z_Eta=cell(1,k);
    for i=1:k
    eta=data.eta_L+(data.eta_U-data.eta_L)*(i-1)/9800;
    [invR_Z_Eta{i},logdetR_Z_Eta(i),R_Z_Eta{i}]=conditional_eta(eta,data);
    end
    data.invR_Z_Eta=invR_Z_Eta;
    data.logdetR_Z_Eta=logdetR_Z_Eta;
    data.R_Z_Eta=R_Z_Eta;
elseif basis_type==5
    invR_Z_Eta=cell(99,99);
    logdetR_Z_Eta=zeros(99,99);
    R_Z_Eta=cell(99,99);
    for i=1:99
        for j=1:99     
            eta=[data.eta_L+(data.eta_U-data.eta_L)*(i-1)/98,data.eta_L+(data.eta_U-data.eta_L)*(j-1)/98];
            [invR_Z_Eta{i,j},logdetR_Z_Eta(i,j),R_Z_Eta{i,j}]=conditional_eta(eta,data);
        end
    end
    data.invR_Z_Eta=invR_Z_Eta;
    data.logdetR_Z_Eta=logdetR_Z_Eta;
    data.R_Z_Eta=R_Z_Eta;
elseif basis_type==7
    k=9801;
    invR_Z_Eta=cell(1,k);
    logdetR_Z_Eta=zeros(1,k);
    R_Z_Eta=cell(1,k);
    for i=1:k
    eta=data.eta_L+(data.eta_U-data.eta_L)*(i-1)/9800;
    [invR_Z_Eta{i},logdetR_Z_Eta(i),R_Z_Eta{i}]=conditional_eta(eta,data);
    end
    data.invR_Z_Eta=invR_Z_Eta;
    data.logdetR_Z_Eta=logdetR_Z_Eta;
    data.R_Z_Eta=R_Z_Eta;
else 
    data.invR_Z_Eta=eye(M);
    data.logdetR_Z_Eta=0;
    data.R_Z_Eta=eye(M);
end
if data.isotropic==0
    ind=data.ind;
    Inv_Sigma_eta=[];
    for i=1:201
        invs2=data.invR_Z_Eta{ind(i,1),ind(i,2)};
        logs2=data.logdetR_Z_Eta(ind(i,1),ind(i,2));
        Inv_Sigma_eta=[Inv_Sigma_eta,invs2];
        Log_Det_Sigma_eta(i)=logs2;
    end
        data.Inv_Sigma_eta=Inv_Sigma_eta;
    data.Log_Det_Sigma_eta=Log_Det_Sigma_eta;
else
    if data.fix_eta==0
        ind=[1:100:9800,9801];
        Inv_Sigma_eta=[];
        for i=1:length(ind)
            invs2=data.invR_Z_Eta{ind(i)};
            logs2=data.logdetR_Z_Eta(ind(i));
            Inv_Sigma_eta=[Inv_Sigma_eta,invs2];
            Log_Det_Sigma_eta(i)=logs2;
        end
            data.Inv_Sigma_eta=Inv_Sigma_eta;
    data.Log_Det_Sigma_eta=Log_Det_Sigma_eta;
    end
end
end
function [invR,logdetR,Sigma_z]=conditional_eta(eta,data)
R=data.R;
Sig=R(data.X,data.Y,eta);
Sigma_z=data.KL_trans'*Sig*data.KL_trans;
%[S,V,D]=svd(Sigma_z);
[invR,logdetR]=invandlogdet(Sigma_z);
end


function data=KLbasislower(data,M)
eta=data.eta_L;
if data.isotropic==0
   eta=[eta,eta];
end
R_L=data.R(data.X,data.Y,eta);
[s,v,d]=svd(R_L);
KL=s(:,1:M);
data.eta=eta;
data.KL=KL*(v(1:M,1:M)^(1/2));
KL_full=s*(v.^0.5);
data.KL_full=KL_full;
data.eig=diag(v);
end

function data=KLbasiscustmize(data,M,eta)

if data.isotropic==0
    if length(eta)==1
        eta=[eta,eta];
    end
else
    if length(eta)==2
        eta=eta(1);
    end
end
R_L=data.R(data.X,data.Y,eta);

[s,v,d]=svd(R_L);
KL=s(:,1:M);
data.eta=eta;
data.KL=KL*(v(1:M,1:M)^(1/2));
KL_full=s*(v.^0.5);
data.KL_full=KL_full;
data.eig=diag(v);
end

function data=KLbasisupper(data,M)
eta=data.eta_U;
if data.isotropic==0
   eta=[eta,eta];
end
R_L=data.R(data.X,data.Y,eta);
[s,v,d]=svd(R_L);
KL=s(:,1:M);
data.eta=eta;
KL_full=s*(v.^0.5);
data.KL_full=KL_full;
data.KL=KL*(v(1:M,1:M)^(1/2));
data.eig=diag(v);
end

function data=KLbasiscentral(data,M)
eta=(data.eta_U+data.eta_L)/2;
if data.isotropic==0
   eta=[eta,eta];
end
R_L=data.R(data.X,data.Y,eta);
[s,v,d]=svd(R_L);
KL=s(:,1:M);
KL_full=s*(v.^0.5);
data.KL_full=KL_full;
data.eta=eta;
data.KL=KL*(v(1:M,1:M)^(1/2));
data.eig=diag(v);
end
