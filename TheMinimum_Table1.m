Mc5v=zeros(1,5);
r0l=zeros(1,5);
r0mv=zeros(1,5);
r0u=zeros(1,5);

zbbhv=[0.0213,0.3056,0.0172,0.0033,0.0422];
m1v=[12,183,1.51,0.44,1.46];
m2v=[7,1.5,1.26,0.12,0.04];
pbbhv=[1.05,12.08,14.1,15.9,1.2];
e0v=[0,0.7,0.13,0.42,0];
% To reproduce Table 1, last four columns
name={'3C 66B','OJ 287','NGC 5548','NGC 4151','Mrk 231'};
for k=1:5
    m1=m1v(k);
    m2=m2v(k);
    zbbh=zbbhv(k);
    pbbh=pbbhv(k);
    e0=e0v(k);
    [Mc1,Tc1,r0n,Anaive]=AyrNaivefun(m1,m2,zbbh,pbbh,e0);
    [Mc,Tc,rate,r0m,Pm15,Ayrm,AyrL]=AyrLowerfun(m1,m2,zbbh,pbbh,e0);
    fprintf('Name     Mc(1e8 Msun)  Tc (Myr)    r0naive  Ayrnaive\n');
    fprintf('%s   %4.2e      %4.2e       %4.2e       %4.2e \n',name{k},Mc1,1e3*Tc1,r0n,Anaive);
    fprintf('----------------------------\n');
    
    fprintf('Name     Mc(1e8 Msun)  Tc (Myr)    r0 (Mpc^-3Gyr^-1) Ayr (1e-16) AyrLow (1e-16)\n');
    fprintf('%s   %4.2e      %4.2e       %4.2e       %4.2f       %4.2f \n',name{k},Mc,1e3*Tc,r0m(2),1e16*Ayrm,1e16*AyrL);
    fprintf('----------------------------\n');
    Mc5v(k)=1e8*Mc;
    r0l(k)=r0m(1);
    r0mv(k)=r0m(2);
    r0u(k)=r0m(3);
end
%% To reproduce Fig. 2
% Mc5v=1e8.*[7.92,10.2312,1.2,0.19,0.17];
sigmc=1e8.*[3.7,0.43,0.37,0,0];
% r0l=[0.0117,2.27e-7,9.9e-7,8.7e-6,1.8e-6];
% r0mv=[0.102,1.976e-6,8.58e-6,7.5e-5,1.59e-5];
% r0u=[0.441,8.547e-6,3.71e-5,3.25e-4,6.89e-5];
name={'3C 66B','OJ 287','NGC 5548','NGC 4151','Mrk 231'};
xdx=1.2.*Mc5v;
figure
xL=1e7;
xU=1e10;
yL=2.5e-5; % galaxy merger rate, lower
yU=6.3e-4; % galaxy merger rate, upper
% see Conselice 2014 (https://arxiv.org/abs/1403.2783)

fillhandle=fill([xL,xU,xU,xL],[yL,yL,yU,yU],'blue');
set(fillhandle,'EdgeColor','k','FaceAlpha',0.45,'EdgeAlpha',0.45);
hold on
%loglog(Mc,ratdat(:,1),'b*')
% errorbar(x,y,yneg,ypos,xneg,xpos)
errorbar(Mc5v,r0mv,r0mv-r0l,r0u-r0mv,sigmc,sigmc,'o','MarkerSize',15,'linewidth',3)
%legend('NGC 5548','Mrk 231','OJ 287','PG 1302','NGC 4151','Ark 120')
text(xdx, r0mv, name,'FontName','Times New Roman','FontSize',20);
xlabel('$M_c$ $(M_{\odot})$'), ylabel('$r_0$ ($\rm{Mpc}^{-3}Gyr^{-1}$)'), axis xy;
set(gca,'XLim',[1e7 1e10], 'xscale', 'log', 'yscale', 'log')
grid on
%%
