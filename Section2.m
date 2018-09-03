%% reproduce numbers in Section 2
G = 6.67408.*10.^(-11);
C = 299792458;
solar_mass = 1.98855e30; 
h=67.8; % Planck 2016 cosmology (arXiv:1502.01589)
Omegam=0.308;
OmegaL=0.692; 
Mpc=1e6 * 3.08567758149137e16 ;
yr=365.25*86400;
Gyr=1e9*yr;
H0 = h*1000/Mpc;
f0=1/yr;
etam=0.25; % maximum symmetric mass ratio when m1=m2
Ccst0=(Mpc^(-3))*(4/3)*(pi^(-1/3))*(C^(-2))*((G*solar_mass)^(5/3))*(f0^(-4/3));
Aconst=sqrt(Ccst0*etam); % constant in Eq. 10

% reference values for Eqs 6-7
r0=1e-4; % (local merger rate density)
N0=1e-4; % (local number density of SMBHs)
Mc1=1e8;

zcut=2; % cutoff redshift
Ez = @(z) sqrt(Omegam*(1+z).^3+OmegaL);
Ind=0;% index m in e(z)=(1+z)^m
phidr = @(z) (1+z).^(Ind-(4/3))./Ez(z);
phir = integral(phidr,0,zcut);
Ez1 = @(z) (1+z).^(Ind-1)./Ez(z);
phirz = integral(Ez1,0,zcut);
zm34=sqrt(phir/0.63); % factor in Eq 6
% 0.8265 for m=-1, 1.2638 m=1
zm13=sqrt(phir/(0.86*phirz)); % Eq 7
% 1.01 for m=-1, 0.987 m=1
Im43=phir;
Im431=phir/phirz;

hcsq=(Mc1^(5/3))*(r0/H0)*Gyr^(-1);
% 4.8209e-16, using r0
Aref6=sqrt(Ccst0*hcsq*Im43); % Eq 6
hcsq=(Mc1^(5/3))*N0;
% 1.4869e-16, using N0
Aref7=sqrt(Ccst0*hcsq*Im431); % Eq 7

% numbers in Phinney 2001 for SMBBH background, page 4
Mc=(10^(7.8))*2*(0.25^(3/5));
ArefP1=3e-24*(Mc^(5/6))*((f0/1e-3)^(-2/3))*sqrt(r0);
%
fprintf('Eq6:Constant    Ialpha \n');
fprintf('%4.2e         %4.2f\n',Aref6,Im43);
fprintf('----------------------------\n');
fprintf('Eq7:Constant    Ialpha \n');
fprintf('%4.2e         %4.2f\n',Aref7,Im431);
fprintf('----------------------------\n');
%%