%% Section 3
% load files for black hole mass functions; see Appendix and references therein for details
fname={'Marconi04.txt','Li2011.txt','Shankar13_red_line.txt','Ueda2014.txt','Mutlu2016.txt','Sesana13.txt'};
Legendname={'Marconi et al 2004','Li et al. 2011','Shankar 2013','Ueda et al. 2014','Mutlu-Pakdil et al. 2016','Sesana 2013'};
ColV={'b','r','g','m','k','c'};
Nmodel=6;
% number of Monte-Carlo simulations to obtain the probability distribution of Amax
Nsim=10000;
DatNMc0=zeros(Nmodel,4);
Ayr95=zeros(3,Nmodel);

% factors in Eq (10)
Aref=1.7228e-21; % constant
zcut=2; % redshift integration limit
Ez1 = @(z) (1+z).^(-1/3);
zm13=(1.1^(5/3))*integral(Ez1,0,zcut);
%
AmaxV=linspace(log10(2e-16),log10(6e-15),128);
Aymv=10.^(AmaxV);
pdfv=zeros(size(Aymv));
tic
for k=1:Nmodel
    file_name=fname{k};
    BHMFmut=load(file_name);
    if k==2
        m=BHMFmut(:,1);
        mf=BHMFmut(:,2);
        mfU=mf+0.2; % 20% in log(BHMF) uncertainty
        mfL=mf-0.2;    
    else
        m=BHMFmut(:,1);
        mf=BHMFmut(:,2);
        mfU=BHMFmut(:,3);
        mfL=BHMFmut(:,4);
    end
    mmax=max(m);
    m0=7:0.01:mmax;
    mv=10.^(m0);
    mf1=interp1(m,mf,m0,'pchip');
    mfU1=interp1(m,mfU,m0,'pchip');
    mfL1=interp1(m,mfL,m0,'pchip');
    mfe=(mfU1-mfL1)./2;
    
    N0=trapz(m0,10.^(mf1));
    Md0=trapz(m0,(10.^m0).*(10.^mf1));
    M0=Md0/N0;
    Mca=M0*(0.25^(3/5)); % 0.2096 is the mean for q in [0.1,1]
    DatNMc0(k,:)=[1e3*N0,1e-5*Md0,1e-8*M0,1e-7*Mca];
    js1=trapz(m0,(10.^(mf1)).*(mv.^(5/3)));
    Ayr95(1,k)=sqrt(zm13)*Aref.*sqrt(js1);
    % Ayr95(1,k) mean Amax using mean BHMF
    % Ayr95(2,k) median Amax
    % Ayr95(3,k) 1-sigma lower Amax
    % Ayr95(4,k) 1-sigma upper Amax
    [AyrV,pdf,Ayr95(2,k),Ayr95(3,k),Ayr95(4,k)]=CompAmaxfun(m0,mf1,mfe,Nsim);
    if k==1
        figure;semilogx(AyrV,pdf,ColV{k},'linewidth',2)
        %set(gca,'Xlim',[1 12])
    else
        hold on;semilogx(AyrV,pdf,ColV{k},'linewidth',2)
        %set(gca,'Xlim',[1 12])
    end
    pdfv=pdfv+interp1(AyrV,pdf,Aymv,'spline',0);
end
%
pdfnm=pdfv./trapz(AmaxV,pdfv);
cdfnm=cumtrapz(AmaxV,pdfnm)+linspace(0, 1, length(pdfnm))*1e-10;
hold off
%legend('Marconi04','Shankar04','Shankar13','Ueda14','Mutlu-Pakdil16','Sesana13')
legend(Legendname{1:Nmodel})
set(gca,'Xlim',[1e-16 1e-14])
grid on
% Fig. 1 in the paper + Sesana 2013 Amax
%
AyrmaxUu=10^(interp1(cdfnm,AmaxV,0.95));
% 2.4e-15, upper limit assuming 5 BHMF models equally probable

%% convert Amax to A_optimistic
BHMFz0=load('U14_z0.txt');
% this is somewhat different from Ueda2014.txt
BHMFz1=load('U14_z1.txt');
BHMFz2=load('U14_z2.txt');
zcut=2;
Ez1 = @(z) (1+z).^(-1/3);
phirz = integral(Ez1,0,zcut);

m0=(7:0.01:10)';
z=0:zcut;
z0=0:0.1:zcut;
Phim0=interp1(BHMFz0(:,1),BHMFz0(:,2),m0,'pchip');
Phim1=interp1(BHMFz1(:,1),BHMFz1(:,2),m0,'pchip');
Phim2=interp1(BHMFz2(:,1),BHMFz2(:,2),m0,'pchip');
PhiMat1=[Phim0,Phim1,Phim2];

m=10.^m0;
m1=m.^(5/3);
mf1=10.^(Phim0);
jis1=trapz(m0,mf1.*m1);
Ayr2=sqrt(phirz*jis1); % numerator in Eq. 11

[X,Y] = meshgrid(z,m0);
[X0,Y0] = meshgrid(z0,m0);
PhiMat0=interp2(X,Y,PhiMat1,X0,Y0,'cubic');

zm13=(1+z0).^(-1/3);% 
zjs=zeros(size(z0));
for k=1:length(z0)
    mf1=PhiMat0(:,k);
    zjs(k)=trapz(m0,10.^(mf1).*m1);
end
js1=trapz(z0,zm13.*zjs);

AyrUedaz=sqrt(js1); % denominator in Eq. 11
%
Eq = @(q) q./(1+q).^2; % definition of symmetric mass ratio, q=m2/m1
etamn = integral(Eq,0,1);
% mean symmetric mass ratio for uniform q in [0,1]

Az=Ayr2/AyrUedaz; % =1.1, Eq.11
feta=sqrt(0.25/etamn); % =1.14
frad=(1.1/1.05)^(5/6); % =1.04
Max2opt=1/(Az*feta*frad); % =0.77

Ayroptim=Max2opt.*Ayr95(3:4,6);
% to compare with Sesana 2013 fiducial model

