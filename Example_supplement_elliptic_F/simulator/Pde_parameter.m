function [True,data]=Pde_parameter(eta_true,z_true,sigma_e,isotropic)
% check
% this function is for generating the true value, observations, determine
% the PDE geometry and parameters.
% sigma_e the propotion of error added into the observations.
% all global variables are PDE parameters.
% Warning: all units of length in this simulation is meter (instead of km). e.g. 
% the unit for x, y are m;
% the unit for T is m^2/day;
% the unit for h(x) is m/day;
% the unit for PDE solution u is m;
% the unit for Dirichlet BC is m;
% ...
global T x_max y_max r hmax
data.eta_L=0.01;
data.eta_U=0.99;
if nargin>3
    data.isotropic=isotropic;%0 means input GP is anisotropic, 1 means input GP is isotropic.
else
    data.isotropic=1;
end
x_max=1000;
y_max=1000;
T=exp(3.8);% transmissivity, the unit is m^2/day.
r=60;  
hmax=50; % the largest mesh length allowed in generating mesh
[x_obs,y_obs]=meshgrid(0.125:0.125:0.875); %observation location, i.e. position of observation wells
x_obs=x_obs*x_max;
y_obs=y_obs*y_max;
pdemodel=settingforpde(hmax);
data.D_in=[x_obs(:),y_obs(:)];
%% Parameter setting
if data.isotropic==1
    True.eta_true=eta_true;
    R=@(X,Y,eta) exp(log(eta)*(((X(:)-X(:)').^2+(Y(:)-Y(:)').^2)/(x_max*y_max)));
else  %anisotropic
    if length(eta_true)==1
        True.eta_true=[eta_true,eta_true];
    else
        True.eta_true=eta_true;
    end

    R=@(X,Y,eta) exp((log(eta(1))*((X(:)-X(:)').^2)+log(eta(2))*((Y(:)-Y(:)').^2))/(x_max*y_max));
end
True.z_true=z_true;
True.pdemodel=pdemodel;
True.sigma_e=sigma_e;
True.x_obs=x_obs(:)';
True.y_obs=y_obs(:)';
True.R=R;
True.hmax=hmax;
True=PRE_SETUP(True);
%data.R_int=R_int;
data.obs=True.obs;
data.X=True.X;
data.Y=True.Y;
data.sigma_e=sigma_e;
data.x_obs=x_obs(:)';
data.y_obs=y_obs(:)';
data.R=R;
data.hmax=hmax;
data.explain_var=cumsum(True.eig)/sum(True.eig);
end

function True=PRE_SETUP(parameter_setting)
    True=parameter_setting;
    True.X=True.pdemodel.Mesh.Nodes(1,:);
    True.Y=True.pdemodel.Mesh.Nodes(2,:);
    Sigma=True.R(True.X,True.Y,True.eta_true);
    [S,V,~]=svd(Sigma);
    M_true=find(cumsum(diag(V))./sum(diag(V))>0.95, 1);
    N=length(True.Y);
    
    %V=eye(size(V,1))*sum(sum(V))/M_true;%%%%%%%%%%%%%%%%%%%%%%%
    
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
function y=settingforpde(hmax)
global x_max y_max
%%
    model = createpde();
    R1 = [3,4,0,0,x_max,x_max,y_max,0,0,y_max]';
    
%     [x_well_location, y_well_location]=meshgrid([0.25:0.25:0.75]);
%     x_well_location=x_well_location*x_max;
%     y_well_location=y_well_location*y_max;
%     x_well_location=x_well_location(:);
%     y_well_location=y_well_location(:);
%     num_well=length(x_well_location);
     geom = [R1];
     sf = 'R1';
%     for i=1:num_well
%         x_well1_min=x_well_location(i)-0.5;
%         x_well1_max=x_well_location(i)+0.5;
%         y_well1_min=y_well_location(i)-0.5;
%         y_well1_max=y_well_location(i)+0.5;
%         R2 = [3,4,x_well1_min,x_well1_min,x_well1_max,x_well1_max,y_well1_max,y_well1_min,y_well1_min,y_well1_max]';
%         R3 = [3,4,x_well1_min-0.01,x_well1_min-0.01,x_well1_max+0.01,x_well1_max+0.01,y_well1_max+0.01,y_well1_min-0.01,y_well1_min-0.01,y_well1_max+0.01]';
%         geom=[geom,R2,R3];
%         sf=[sf,'+R',num2str(2*i),'+R',num2str(2*i+1)];
%     end
    % Names for the two geometric objects
    ns = (char('R1'))';
    % Set formula
    %sf = 'R1+R2+R3';
    gd= decsg(geom,sf,ns);
    geometryFromEdges(model,gd);  
    generateMesh(model,'Hmax',hmax);
    y=model;
    %pdegplot(model,'EdgeLabels','on','FaceLabels','on')
end
