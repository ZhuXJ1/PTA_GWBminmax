function [Mc,Tc,rate,Anaive]=AyrNaivefun(m1,m2,zbbh,pbbh,e0)
% Input: SMBBH parameters
% m1 and m2 (component masses, in 1e8 Msun)
% zbbh (redshift)
% pbbh (observed orbital period, in yr)
% e0 (orbital eccentricity)

% Output: Mc - chirp mass (in 1e8 Msun),
% Tc - coalescence time (in Gyr)
% rate - naive merger rate estimate (see Eq. 12 in the paper, in Mpc^-3Gyr^-1)
% Anaive - naive estimator of Ayr, insert rate and Mc into Eq. 6 in the paper

% Note: this is the code to use when you have a new SMBBH candidate and
% want to know its implication for the GWB. rate and Anaive will tell you
% if it would imply a creazily low or high merger rate and GWB amplitude.

C = 299792458;
solar_mass = 1.98855e30; 
G = 6.67408e-11;
yr=365.25*86400;
Gyr=1e9*yr;
Omegam=0.308;
OmegaL=0.692;
cc=299792.458;
h=67.8;
zv=0:0.001:1;

Ez = sqrt(Omegam.*(1+zv).^3+OmegaL);
phidr = 1./Ez;
drz=(cc/h).*cumtrapz(zv,phidr);
dVdz = 4*pi*(cc/h).*phidr.*(drz.^2);
dVz = cumtrapz(zv,dVdz); % Mpc^3, comoving volume

dedt = @(e) e.^(29/19).*(1+(121.*(e.^2)./304)).^(1181/2299).*((1-e.^2).^(-3/2));

mbbh=m1+m2;
qbbh=m2./m1;
yita=qbbh./((1+qbbh).^2);
Mc=mbbh.*yita.^(3/5);
M_c=1e8*solar_mass.*Mc;
w0=2*pi*(1+zbbh)./(yr.*pbbh);
if e0==0
    Tc=(5/256)*(C^5).*((G.*M_c).^(-5/3)).*(w0.^(-8/3))/Gyr;
else
    chi0=(1-e0.^2).*(e0.^(-12/19)).*((1+121*(e0.^2)/304).^(-870/2299));
    js1=integral(dedt,0,e0);
    Tc=(15/304)*(C^5).*((G.*M_c).^(-5/3)).*(w0.^(-8/3)).*(chi0.^4).*js1/Gyr;
end
tj=interp1(zv,dVz,zbbh,'spline'); % comoving volume, which is actually equivalent to (4\pi d^3)/3
rate=1./tj./Tc; % Eq. 12 in the paper

hc0=4.8e-16.*((Mc).^(5/6));
% Eq. 6 in the paper, 4.8e-16 from Section2.m
Anaive=hc0.*((1e4.*rate).^(1/2));
return
