function y=Simulator(Z,data)
global T time time_obs S 
%computer experiment, simulator, 
% Warning: all units involves length in this simulation is meter (instead of km). e.g. 
% the unit for x, y are m;
% the unit for T is m^2/s;
% the unit for h(x) is m/s;
% the unit for PDE solution u is m;
% the unit for Dirichlet BC is m;
% ....
%input is the calibration parameter, in this problem stands for the functional input
%data contains experimental setting
X=data.X; % unit m 
Y=data.Y; % unit m
input=exp(Z*data.KL'); % unit: m^2/s
model=settingforpde(data.hmax);%construct the PDE geometry
applyBoundaryCondition(model,'dirichlet','Edge',[1,3],'u',45);
applyBoundaryCondition(model,'neumann','Edge',[2,4],'g',0);
C= scatteredInterpolant(X',Y',input');
K=@(location, state) C(location.x,location.y)*T;%this is the setup of coefficient "c", i.e. Transmissivity.
% the function form, m u_tt + d u_t-\nabla(c\nabla u)+a u=f
specifyCoefficients(model,'m',0,'d',S,'c',K,'a',0,'f',@fsource);
setInitialConditions(model,45);
result = solvepde(model,time);
y=interpolateSolution(result,data.x_obs(:),data.y_obs(:),time_obs);
y=y(:)';
end
function y=fsource(location, state)%this is the setup of source term
global r x_max y_max
    [x_well_location, y_well_location]=meshgrid([0.25:0.25:0.75]);
    x_well_location=x_well_location*x_max;
    y_well_location=y_well_location*y_max;
    x_well_location=x_well_location(:);
    y_well_location=y_well_location(:);
    num_well=length(x_well_location);
    y=0;
    for i=1:num_well
        D=((location.x-x_well_location(i)).^2+(location.y-y_well_location(i)).^2).^(1/2);
        y=y-r*normpdf(D/10)/10^2;
    end
end
function y=settingforpde(hmax)
global x_max y_max
%%
    model = createpde();
    R1 = [3,4,0,0,x_max,x_max,y_max,0,0,y_max]';
    geom = [R1];
    sf = 'R1';
    ns = (char('R1'))';
    gd= decsg(geom,sf,ns);
    geometryFromEdges(model,gd);  
    generateMesh(model,'Hmax',hmax);
    y=model;
end
