eta=[0.01*ones(1,10),0.5*ones(1,10),0.99*ones(1,10)];
ETA=repmat(eta,1,10);
%% plot the relative error of differen priro specification.
load('ESS_prior.mat');
test_res(ESS_post,ETA)
ERR_prior=ESS_post;
Hgcf=figure('color','w');
l1={'Unif-GP','L-GP','M-GP','U-GP'};
%ERR_prior=ESS_post(:,[1,2,3,4]);
%ERR_prior(:,1)=min(ESS_post(:,1),ESS_post(:,5));
subplot(1,3,1)
l=l1;
%subplot(1,2,1)
boxplot(ERR_prior(ETA==0.01,1:4),'Labels',l )
set(gca,'yTick',[0:0.25:1.5]);
ylabel('Relative $L_2$ error','Interpreter','latex','Fontname', 'Times New Roman','FontSize',12)
ylim([0 1.5])
title('(a) True value of $\eta =0.01$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',12)
subplot(1,3,2)
%subplot(1,2,1)
boxplot(ERR_prior(ETA==0.5,1:4),'Labels',l )
ylabel('Relative $L_2$ error','Interpreter','latex','Fontname', 'Times New Roman','FontSize',12)
ylim([0 .75])
title('(b) True value of $\eta =0.5','Interpreter','latex','Fontname', 'Times New Roman','FontSize',12)
subplot(1,3,3)
%subplot(1,2,1)
boxplot(ERR_prior(ETA==0.99,1:4),'Labels',l )
ylabel('Relative $L_2$ error','Interpreter','latex','Fontname', 'Times New Roman','FontSize',12)
ylim([0 .3])
title('(c) True value of $\eta=0.99 $','Interpreter','latex','Fontname', 'Times New Roman','FontSize',12)
set(Hgcf,'position',[0 0 1200 240]) 
print(Hgcf,'-dtiff','-r660','FigureF2.emf');
%% boxplot of different design size
load('ESS_criteria.mat');
ERR_criteria=ESS_post;
ERR_criteria_compare=[ERR_criteria(:,1)-ERR_criteria(:,2),ERR_criteria(:,1)-ERR_criteria(:,3),ERR_criteria(:,2)-ERR_criteria(:,3)];
Hgcf=figure('color','w');
l={'WPV-VL','WPV-AEI','VL-AEI'};
subplot(1,2,1)
boxplot(ERR_criteria_compare,'Labels',l); hold on
plot([0,4],[0,0],'k--')
set(gca,'yTick',[-0.5:0.25:0.5]);
ylabel('Difference in relative $L_2$ error','Interpreter','latex','Fontname', 'Times New Roman','FontSize',12)
%xlabel('Design criterion')
ylim([-0.5 .5])
title(' (a) ','Interpreter','latex','Fontname', 'Times New Roman','FontSize',12)
load('ESS_size.mat');
ERR_size=ESS_post;
l2={'0','20','50'};
subplot(1,2,2)
boxplot(ERR_size,'Labels',l2)
ylabel('Relative $L_2$ error','Interpreter','latex','Fontname', 'Times New Roman','FontSize',12)
xlabel('Size of follow-up design','Interpreter','latex','Fontname', 'Times New Roman','FontSize',12)
ylim([0 .75])
title('(b)','Interpreter','latex','Fontname', 'Times New Roman','FontSize',12)
set(Hgcf,'position',[0 0 900 300]) 
print(Hgcf,'-dtiff','-r660','FigureF1.emf');
%% different truncation number
load('ESS_M.mat');
ESS_M=ESS_post;
Hgcf=figure('color','w');
label={'3','9','12','15'};
boxplot(ESS_M,'Labels',label)
ylim([0,1])
ylabel('Relative $L_2$ error','Interpreter','latex','Fontname', 'Times New Roman','FontSize',12)
%title('Compare of relative $L_2$ error with different truncation number $M$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',12)
xlabel('$M$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',12)
set(Hgcf,'position',[0 0 400 200]) 
print(Hgcf,'-dtiff','-r660','FigureF3.emf');
m39=signrank(ESS_M(:,1),ESS_M(:,2));
m312=signrank(ESS_M(:,1),ESS_M(:,3));
m315=signrank(ESS_M(:,1),ESS_M(:,4));
%%  analysis of MCMC
load('RESULT_MCMC.mat');
fs=13;
RESULT={RES_opt{1}{2},RES_opt{2}{2},RES_opt{3}{2},RES_opt{4}{2}};
chain=RESULT{1}.chain;
[acf1,lags1,bounds1,h1]=autocorr(chain(:,1),50);
close(gcf);
[acf2,lags2,bounds2,h2]=autocorr(chain(:,10),50);
close(gcf);
Hgcf=figure('color','w');
set(Hgcf,'position',[0 0 1250 220])
subplot(1,3,1)
auto_custom(acf1,acf2,lags1,fs)
subplot(1,3,2)
ksdensity(RESULT{1}.chain(10001:50:end,10),'Support',[0.01,0.99])
hold on 
%plot([0.5,0.5],[0,2],'--','LineWidth',2)
%legend('density estimate','true value','Location','northwest')
%legend('boxoff');
xlim([0 1.01])
xlabel('$\eta$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',fs)
ylabel('Density estimate','Interpreter','latex','Fontname', 'Times New Roman','FontSize',fs)
title('(b)','Interpreter','latex','Fontname', 'Times New Roman','FontSize',fs);
label={'Unif-GP','L-GP','M-GP','U-GP'};%['\eta=' num2str(data.eta)]
set(gca,'Fontname', 'Times New Roman','FontSize',fs)
subplot(1,3,3)
ERR_l2=[RESULT{1}.Err2;RESULT{2}.Err2;RESULT{3}.Err2;RESULT{4}.Err2]';
boxplot(ERR_l2(10001:50:end,:)/(mean(RESULT{1}.True.f_input.^2))^(1/2),label)
set(gca,'Fontname', 'Times New Roman','FontSize',fs,'XTickLabelRotation',30)
ylabel(' Relative $L_2$ error','Interpreter','latex','Fontname', 'Times New Roman','FontSize',fs)
title('(c)','Interpreter','latex','Fontname', 'Times New Roman','FontSize',fs);
set(gca,'position',[0.72 0.26 0.2 0.6]) 
set(Hgcf,'position',[0 0 1050 220]) 
print(Hgcf,'-dtiff','-r660','FigureF4.emf');

%% sample of functional input at diagonal line 
Visualization_of_Result_mean_l(RESULT)
%% effective sample size
mESS(1) = MultiMCMCESS(RESULT{1}.chain(10001:end,:));
mESS(2) = MultiMCMCESS(RESULT{2}.chain(10001:end,[1:9,11]));
mESS(3) = MultiMCMCESS(RESULT{3}.chain(10001:end,[1:9,11]));
mESS(4) = MultiMCMCESS(RESULT{4}.chain(10001:end,[1:9,11]));
%%
function y=MultiMCMCESS(X)
% the batch size is b=floor(n^(1/2)), if n is not a integer multiple of b,
% then ignore the rest, only use n'=b*floor(n/b) sizes to compute the ESS.
Lambda=cov(X);
n=size(X,1);
p=size(X,2);
b=floor(n^(1/2));
a=floor(n/b);
for i=1:a
   Xb(i,:)=mean(X(b*(i-1)+1:b*i,:)); 
end
X_m=mean(X);
X1=(Xb-X_m);
Sigma=b*cov(X1);
det(Lambda)
det(Sigma)
y=n*(det(Lambda)/det(Sigma))^(1/p);
end
%%
function y=Visualization_of_Result_mean_l(RESULT)
% RESULT contains 4 elements, 
N=size(RESULT{1}.chain,1);
data=RESULT{4}.data;
global x_max y_max
x_max=1000;
y_max=1000;
True=RESULT{1}.True;
y_L=min(True.f_input)-0.5;
y_U=max(True.f_input)+0.5;
y_tic=[floor(y_L):0.5:ceil(y_U)];
x_d=(0:0.02:1)*x_max;
y_d=(0:0.02:1)*y_max;
Hgcf=figure('color','w');
subplot(1,4,1)
[F,c0]=plot_diagonal(RESULT{1}.data,True,RESULT{1}.chain);
c_MAP=RESULT{1}.Z_MAP*RESULT{1}.data.KL';
C= scatteredInterpolant(data.X',data.Y',c_MAP');
c_MAP=C(x_d,y_d);%this is the setup of coefficient "c"
hold on 
plot(x_d/x_max,F','color',[0.8,0.8,0.8]);
ylim([y_L y_U])
pic1=plot(x_d/x_max,c0,'r--','LineWidth',2);
pic2=plot(x_d/x_max,mean(F),'b-*','MarkerIndices',1:10:length(c_MAP));
set(gca,'ytick',y_tic);

h1=xlabel('$x_1$');
h2=ylabel('$f(\mathbf{x})$');
set(h1,'Interpreter','latex');
set(h2,'Interpreter','latex');
title('(a) Unif-GP','Interpreter','latex','Fontname', 'Times New Roman','FontSize',12)
set(gca,'Fontname', 'Times New Roman','FontSize',12)
AX2(1,1)=gca;

subplot(1,4,2)
[F,c0]=plot_diagonal(RESULT{2}.data,True,RESULT{2}.chain);
c_MAP=RESULT{2}.Z_MAP*RESULT{2}.data.KL';
C= scatteredInterpolant(data.X',data.Y',c_MAP');
c_MAP=C(x_d,y_d);%this is the setup of coefficient "c"
hold on 
plot(x_d/x_max,F','color',[0.8,0.8,0.8]);
ylim([y_L y_U])
pic1=plot(x_d/x_max,c0,'r--','LineWidth',2);
pic2=plot(x_d/x_max,mean(F),'b-*','MarkerIndices',1:10:length(c_MAP));
set(gca,'ytick',y_tic);

%legend([pic1,pic2],'true','MAP')
h1=xlabel('$x_1$');
h2=ylabel('$f(\mathbf{x})$');
set(h1,'Interpreter','latex');
set(h2,'Interpreter','latex');
title('(b) L-GP','Interpreter','latex','Fontname', 'Times New Roman','FontSize',12)
set(gca,'Fontname', 'Times New Roman','FontSize',12)
AX2(1,2)=gca;

subplot(1,4,3)
[F,c0]=plot_diagonal(RESULT{3}.data,True,RESULT{3}.chain);
c_MAP=RESULT{3}.Z_MAP*RESULT{3}.data.KL';
C= scatteredInterpolant(data.X',data.Y',c_MAP');
c_MAP=C(x_d,y_d);%this is the setup of coefficient "c"
hold on 
plot(x_d/x_max,F','color',[0.8,0.8,0.8]);
ylim([y_L y_U])
pic1=plot(x_d/x_max,c0,'r--','LineWidth',2);
pic2=plot(x_d/x_max,mean(F),'b-*','MarkerIndices',1:10:length(c_MAP));
set(gca,'ytick',y_tic);

%legend([pic1,pic2],'true','MAP')
h1=xlabel('$x_1$');
h2=ylabel('$f(\mathbf{x})$');
set(h1,'Interpreter','latex');
set(h2,'Interpreter','latex');
title('(c) M-GP','Interpreter','latex','Fontname', 'Times New Roman','FontSize',12)
set(gca,'Fontname', 'Times New Roman','FontSize',12)
AX2(1,3)=gca;

subplot(1,4,4)
[F,c0]=plot_diagonal(RESULT{4}.data,True,RESULT{4}.chain);
c_MAP=RESULT{4}.Z_MAP*RESULT{4}.data.KL';
C= scatteredInterpolant(data.X',data.Y',c_MAP');
c_MAP=C(x_d,y_d);%this is the setup of coefficient "c"
hold on 
plot(x_d/x_max,F','color',[0.8,0.8,0.8]);
ylim([y_L y_U])
pic1=plot(x_d/x_max,c0,'r--','LineWidth',2);
pic2=plot(x_d/x_max,mean(F),'b-*','MarkerIndices',1:10:length(c_MAP));
set(gca,'ytick',y_tic);

legend([pic1,pic2],'True functional input','Posterior mean','FontSize',10)
legend('boxoff');
h1=xlabel('$x_1$');
h2=ylabel('$f(\mathbf{x})$');
set(h1,'Interpreter','latex');
set(h2,'Interpreter','latex');
title('(d) U-GP','Interpreter','latex','Fontname', 'Times New Roman','FontSize',12)
set(gca,'Fontname', 'Times New Roman','FontSize',12)
AX2(1,4)=gca;

set(Hgcf,'units','pixel');
set(Hgcf,'position',[0 0 1000,250]) 
set(Hgcf, 'PaperPositionMode', 'auto');
print(Hgcf,'-dtiff','-r660','FigureF5.emf');
Z_mean1=mean(RESULT{1}.chain(N/6:50:end,1:RESULT{1}.data.M));
Z_mean2=mean(RESULT{2}.chain(N/6:50:end,1:RESULT{2}.data.M));
Z_mean3=mean(RESULT{3}.chain(N/6:50:end,1:RESULT{3}.data.M));
Z_mean4=mean(RESULT{4}.chain(N/6:50:end,1:RESULT{4}.data.M));
e1=[mean(abs(Z_mean1*RESULT{1}.data.KL'-True.f_input)),mean(abs(Z_mean2*RESULT{2}.data.KL'-True.f_input)),mean(abs(Z_mean3*RESULT{3}.data.KL'-True.f_input)),mean(abs(Z_mean4*RESULT{4}.data.KL'-True.f_input))];
e2=[mean((Z_mean1*RESULT{1}.data.KL'-True.f_input).^2),mean((Z_mean2*RESULT{2}.data.KL'-True.f_input).^2),mean((Z_mean3*RESULT{3}.data.KL'-True.f_input).^2),mean((Z_mean4*RESULT{4}.data.KL'-True.f_input).^2)];
e3=[max(abs(Z_mean1*RESULT{1}.data.KL'-True.f_input)),max(abs(Z_mean2*RESULT{2}.data.KL'-True.f_input)),max(abs(Z_mean3*RESULT{3}.data.KL'-True.f_input)),max(abs(Z_mean4*RESULT{4}.data.KL'-True.f_input))];
y=[e1;e2.^(1/2);e3];

end

function [y,c]=plot_diagonal(data,True,chain)
global x_max y_max
N=size(chain,1);
X=data.X;
Y=data.Y;
chain=chain(10001:50:end,:);
n=size(chain,1);
x_d=(0:0.02:1)*x_max;
y_d=(0:0.02:1)*y_max;
for i=1:n
   f=chain(i,1:data.M)*data.KL';
   C= scatteredInterpolant(X',Y',f');
   c= C(x_d,y_d);%this is the setup of coefficient "c"
   F(i,:)=c; 
end
C= scatteredInterpolant(X',Y',True.f_input');
c= C(x_d,y_d);%this is the setup of coefficient "c"
%c_MAP=MAP_est(F);
y=F;
end


function y=auto_custom(acf1,acf2,lags1,fs)
hold on
for i=1:length(acf1)/2
p1=plot([lags1(2*i-1),lags1(2*i-1)],[0,acf1(2*i-1)],'-.o','MarkerIndices',2,'Color','r','MarkerSize',7);
p2=plot([lags1(2*i-1),lags1(2*i-1)],[0,acf2(2*i-1)],'-s','MarkerIndices',2,'Color','b','MarkerSize',7);
end
title('(a)','Interpreter','latex','Fontname', 'Times New Roman','FontSize',fs);
xlim([0,lags1(end)]);
xlabel('Lag')
plot([0,50],[0,0],'k')
ylabel('Sample autocorrelation','Fontname', 'Times New Roman','FontSize',fs)
legend([p1,p2],{'Autocorrelation of $\zeta_1$','Autocorrelation of $\eta$'},'Interpreter','latex','Fontname', 'Times New Roman')
legend('boxoff');

ylim([-0.1,1.1])
set(gca,'Fontname', 'Times New Roman','FontSize',fs)
end
