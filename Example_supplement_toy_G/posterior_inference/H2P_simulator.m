function y=H2P_simulator(data,start_set)
%%
% global optimazation to find posterior mode. 
global posterior_opt
posterior_opt=data;
posterior_opt.n=length(data.obs);
M=data.M;
%M=10;
posterior_opt.N=size(data.X,1);
    if size(start_set,2)==M+2
       start_set(:,M+1)=[]; 
    elseif size(start_set,2)==M
        start_set(:,M+1)=rand(size(start_set,1),1)*(data.eta_U-data.eta_L)+data.eta_L;
    end
    if data.fix_eta==0
        if data.isotropic==0
            options_patternsearch = optimoptions('patternsearch','Display','off','MaxIterations',3000*M);
            options = optimoptions('particleswarm','Display','off','MaxIterations',3000*M,'HybridFcn',@patternsearch);
            [x01,fval0,exitflag0]= particleswarm(@ESS,data.M,[-10*ones(1,data.M)],[10*ones(1,data.M)],options);
            [x2,fval2,exitflag2] = patternsearch(@minuslogposterior_marginal_aniso,[x01],[],[],[],[],[-10*ones(1,data.M)],[10*ones(1,data.M)],[],options_patternsearch);
            problem = createOptimProblem('fminunc','x0',zeros(1,data.M),'objective',@minuslogposterior_marginal_aniso);
            tpoints=CustomStartPointSet(start_set(:,1:data.M));
            ms=MultiStart('StartPointsToRun','all');
            [x02,fval0,exitflag0]=run(ms,problem,tpoints);
            [x3,fval3,exitflag3] = patternsearch(@minuslogposterior_marginal_aniso,x02,[],[],[],[],[-10*ones(1,data.M)],[10*ones(1,data.M)],[],options_patternsearch);
            initial_swarm_mat=[x01;x02];
            options = optimoptions('particleswarm','Display','off','SwarmSize',10*(M+1)+2,'InitialSwarmMatrix',initial_swarm_mat,'HybridFcn',@patternsearch);
            [x1,fval1,exitflag1] = particleswarm(@minuslogposterior_marginal_aniso,data.M,[-10*ones(1,data.M)],[10*ones(1,data.M)],options);
        else
            options_patternsearch = optimoptions('patternsearch','Display','off','MaxIterations',3000*M);
            options = optimoptions('particleswarm','Display','off','MaxIterations',3000*M,'HybridFcn',@patternsearch);
            
            [x01,fval0,exitflag0]= particleswarm(@ESS,M,[-10*ones(1,M)],[10*ones(1,M)],options);
            [x2,fval2,exitflag2] = patternsearch(@minuslogposterior_marginal,[x01],[],[],[],[],[-10*ones(1,M)],[10*ones(1,M)],[],options_patternsearch);
            problem = createOptimProblem('fminunc','x0',zeros(1,M),'objective',@minuslogposterior_marginal);
            tpoints=CustomStartPointSet(start_set(:,1:M));
            ms=MultiStart('StartPointsToRun','all');
            [x02,fval0,exitflag0]=run(ms,problem,tpoints);
            [x3,fval3,exitflag3] = patternsearch(@minuslogposterior_marginal,x02,[],[],[],[],[-10*ones(1,M)],[10*ones(1,M)],[],options_patternsearch);
            initial_swarm_mat=[x01;x02];
            options = optimoptions('particleswarm','Display','off','SwarmSize',10*(M+1)+2,'InitialSwarmMatrix',initial_swarm_mat,'HybridFcn',@patternsearch);
            [x1,fval1,exitflag1] = particleswarm(@minuslogposterior_marginal,data.M,[-10*ones(1,M)],[10*ones(1,M)],options);
        end
    else
        options_patternsearch = optimoptions('patternsearch','Display','off','MaxIterations',3000*M);
        options = optimoptions('particleswarm','Display','off','MaxIterations',3000*M,'HybridFcn',@patternsearch);
        [x01,fval0,exitflag0]= particleswarm(@ESS,data.M,[-10*ones(1,data.M)],[10*ones(1,data.M)],options);
        [x2,fval2,exitflag2] = patternsearch(@minuslogposterior_fix_eta,[x01],[],[],[],[],[-10*ones(1,data.M)],[10*ones(1,data.M)],[],options_patternsearch);
        problem = createOptimProblem('fminunc','x0',zeros(1,data.M),'objective',@minuslogposterior_fix_eta);
        tpoints=CustomStartPointSet(start_set(:,1:data.M));
        ms=MultiStart('StartPointsToRun','all');
        [x02,fval0,exitflag0]=run(ms,problem,tpoints);
        [x3,fval3,exitflag3] = patternsearch(@minuslogposterior_fix_eta,x02,[],[],[],[],[-10*ones(1,data.M)],[10*ones(1,data.M)],[],options_patternsearch);
        initial_swarm_mat=[x01;x02];
        options = optimoptions('particleswarm','Display','off','SwarmSize',10*(M+1)+2,'InitialSwarmMatrix',initial_swarm_mat,'HybridFcn',@patternsearch);
        [x1,fval1,exitflag1] = particleswarm(@minuslogposterior_fix_eta,data.M,[-10*ones(1,data.M)],[10*ones(1,data.M)],options);
    end
    y.x_local=[x1,fval1;x2,fval2;x3,fval3];%;x4,fval4'];
    fval=y.x_local(:,end);
    ind=find(fval==min(fval),1);
    y.opt_ind=ind;
    y.x_opt=y.x_local(ind,1:data.M+1);
    y.x_opt_h2p=y.x_local(1,1:data.M+1);
end

%%
function y=minuslogposterior_marginal(para)
global posterior_opt
sigma_e=posterior_opt.sigma_e;
M=posterior_opt.M;
N=size(posterior_opt.KL,2);
zeta=para(1:M);
m=Simulator([zeta,zeros(1,N-M)],posterior_opt);
logs1=sum(log(sigma_e^2));
minusloglik=logs1/2+0.5*sum((posterior_opt.obs-m).^2./(sigma_e^2));
minuslogprior=-log(prior_marginal_z(zeta,posterior_opt));
y=minusloglik+minuslogprior;
end
function prior=prior_marginal_z(z,data)
M=data.M;
N=size(data.KL,2);
temp=[z,zeros(1,N-M)]*(reshape(data.Inv_Sigma_eta'*[z,zeros(1,N-M)]',data.M,99));
prior=mean(exp(-data.Log_Det_Sigma_eta/2).*((temp/2++data.para_prior.b).^(-data.M/2-data.para_prior.a)));
end

%%
function y=minuslogposterior_marginal_aniso(para)
global posterior_opt
sigma_e=posterior_opt.sigma_e;
M=posterior_opt.M;
zeta=para(1:M);
m=Simulator(zeta,posterior_opt);

logs1=sum(log(sigma_e^2));
minusloglik=logs1/2+0.5*sum((posterior_opt.obs-m).^2./(sigma_e^2));
minuslogprior=-log(prior_marginal_z_aniso(zeta,posterior_opt));
y=minusloglik+minuslogprior;
end

function prior=prior_marginal_z_aniso(z,data)
temp=z*(reshape(data.Inv_Sigma_eta'*z',data.M,201));
prior=mean(exp(-data.Log_Det_Sigma_eta/2).*((temp/2++data.para_prior.b).^(-data.M/2-data.para_prior.a)));
end

%%
function y=minuslogposterior_fix_eta(para)
global posterior_opt 

sigma_e=posterior_opt.sigma_e;
M=posterior_opt.M;
zeta=para(1:M);
m=Simulator(zeta,posterior_opt);

logs1=sum(log(sigma_e^2));
minusloglik=logs1/2+0.5*sum((posterior_opt.obs-m).^2./(sigma_e^2));
minuslogprior=(M/2+posterior_opt.para_prior.a)*log(zeta*zeta'/2+posterior_opt.para_prior.b);
y=minusloglik+minuslogprior;
end

function y=ESS(para)
global posterior_opt
zeta=para(1:posterior_opt.M);
m=Simulator(zeta,posterior_opt);
y=sum((posterior_opt.obs-m).^2);
end

