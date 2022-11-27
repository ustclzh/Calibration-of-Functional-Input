load('RESULT_MCMC.mat');
load('ESS_exact.mat');
boxplot(ESS_post,'Labels',{'Proposed method','No-emulator method','Exact method'});
ylabel('Relative $L_2$ error','Interpreter','latex','Fontname', 'Times New Roman','FontSize',14)
ylim([0 1])
set(gca,'Fontname', 'Times New Roman','FontSize',14);
set(gcf,'units','pixel');
set(gcf,'position',[0 0 600,240]) 
set(gcf, 'PaperPositionMode', 'auto');
print(gcf,'-dtiff','-r660','FigureG1.emf');
%%
Visualization_of_Result({RES_opt{1}{2},RES_opt{2}{2},RES_opt{3}{2}})
function Visualization_of_Result(RESULT)
N=size(RESULT{1}.chain,1);
data=RESULT{1}.data;
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
subplot(1,3,1)
[F,c0]=plot_diagonal(RESULT{1}.data,True,RESULT{1}.chain);
c_MAP=RESULT{1}.Z_MAP*RESULT{1}.data.KL';
hold on 
plot(data.X,F','color',[0.8,0.8,0.8]);
ylim([y_L y_U])
pic1=plot(data.X,c0,'r--','LineWidth',2.2);
pic2=plot(data.X,mean(F),'b-*','MarkerIndices',1:10:length(c_MAP));
set(gca,'ytick',y_tic);
h1=xlabel('$x$');
h2=ylabel('$f(x)$');
set(h1,'Interpreter','latex');
set(h2,'Interpreter','latex');
title('(a) Proposed method','Interpreter','latex','Fontname', 'Times New Roman','FontSize',12)
set(gca,'Fontname', 'Times New Roman','FontSize',12)
AX2(1,1)=gca;
subplot(1,3,2)
[F,c0]=plot_diagonal(RESULT{2}.data,True,RESULT{2}.chain);
c_MAP=RESULT{2}.Z_MAP*RESULT{2}.data.KL';
hold on 
plot(data.X,F','color',[0.8,0.8,0.8]);
ylim([y_L y_U])
pic1=plot(data.X,c0,'r--','LineWidth',2.2);
pic2=plot(data.X,mean(F),'b-*','MarkerIndices',1:10:length(c_MAP));
set(gca,'ytick',y_tic);
h1=xlabel('$x$');
h2=ylabel('$f(x)$');
set(h1,'Interpreter','latex');
set(h2,'Interpreter','latex');
title('(b) No-emulator method','Interpreter','latex','Fontname', 'Times New Roman','FontSize',12)
set(gca,'Fontname', 'Times New Roman','FontSize',12)
subplot(1,3,3)
[F,c0]=plot_diagonal(RESULT{3}.data,True,RESULT{3}.chain);
c_MAP=RESULT{3}.Z_MAP*RESULT{3}.data.KL';
hold on 
plot(data.X,F','color',[0.8,0.8,0.8]);
ylim([y_L y_U])
pic1=plot(data.X,c0,'r--','LineWidth',2.2);
pic2=plot(data.X,mean(F),'b-*','MarkerIndices',1:10:length(c_MAP));
set(gca,'ytick',y_tic);
h1=xlabel('$x$');
h2=ylabel('$f(x)$');
set(h1,'Interpreter','latex');
set(h2,'Interpreter','latex');
title('(c) Exact method','Interpreter','latex','Fontname', 'Times New Roman','FontSize',12)
set(gca,'Fontname', 'Times New Roman','FontSize',12)
AX2(1,1)=gca;
legend([pic1,pic2],'True functional input','Posterior mean','FontSize',12)
set(Hgcf,'units','pixel');
set(Hgcf,'position',[0 0 1200,200]) 
set(Hgcf, 'PaperPositionMode', 'auto');
print(Hgcf,'-dtiff','-r660','FigureG2.emf');
end
function [y,c]=plot_diagonal(data,True,chain)
chain=chain(10001:50:end,:);
n=size(chain,1);
for i=1:n
   f=chain(i,1:data.M)*data.KL';
   F(i,:)=f; 
end
c= True.f_input;
y=F;
end