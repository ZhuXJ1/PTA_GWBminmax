function [Mc,Tc,rate,r0m,Pm15,Ayrm,AyrL]=AyrLowerfun(m1,m2,zbbh,pbbh,e0)
% Input: SMBBH parameters
% m1 and m2 (component masses, in 1e8 Msun)
% zbbh (redshift)
% pbbh (observed orbital period, in yr)
% e0 (orbital eccentricity)

% Output: Mc - chirp mass (in 1e8 Msun),
% Tc - coalescence time (in Gyr)
% rate - naive merger rate estimate (see Eq. 12 in the paper)
% r0m - [1-sigma Lower, median, 1-sigma Upper] rate estimates (in Mpc^-3Gyr^-1)
% Pm15 - probability that the implied Ayr is above 1e-15 (the PTA upper limit)
% Ayrm - median estimate of Ayr
% AyrL - 95% lower bound of Ayr

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
dz=interp1(zv,drz,zbbh,'spline'); % comoving distance
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
tj=interp1(zv,dVz,zbbh,'spline');
rate=1./tj./Tc;

rlv=-8:0.01:3;
% the range of priors for r0, has no effect on posterior distribution
r0v=rate.*10.^(rlv);
nmb=1:0.05:100;
% the range of priors for d_max, has no effect on posterior distribution
dzv=dz.*nmb;
m=length(dzv);
n=length(r0v);
pdat=zeros(m,n);

for j=1:m
    for k=1:n
        tj1=(pi*4/3)*(dzv(j)^3);
        pdat(j,k)=r0v(k)*tj1*Tc*exp(-r0v(k)*tj1*Tc)*(dzv(j)^(-3)); % uniform in d_max prior
    end
end

zjs=zeros(size(r0v));
for k=1:length(r0v)
    mf1=pdat(:,k);
    zjs(k)=trapz(dzv,mf1);
end
js1=trapz(log10(r0v),zjs);
pdatn=pdat./js1;
posr0=zeros(size(r0v));
for k=1:length(r0v)
    mf1=pdatn(:,k);
    posr0(k)=trapz(dzv,mf1);
end

r0v1=log10(r0v);
cdf=cumtrapz(r0v1,posr0)+linspace(0, 1, length(posr0))*1e-10;
r0Low=10.^(interp1(cdf,r0v1,0.16));
r0mn=10.^(interp1(cdf,r0v1,0.5));
r0High=10.^(interp1(cdf,r0v1,0.84));
r0m=[r0Low,r0mn,r0High];

hc0=4.8e-16.*((Mc).^(5/6));
hc287=hc0.*((1e4.*r0v).^(1/2));

pdfr=posr0;
cst=trapz(log10(hc287),pdfr);
pdf1=(1/cst).*pdfr;
grd=log10(hc287);
cdf=cumtrapz(grd,pdf1)+linspace(0, 1, length(pdf1))*1e-10;
Ayrm=10.^(interp1(cdf,grd,0.5,'spline')); % median
AyrL=10.^(interp1(cdf,grd,0.05,'spline')); % 95% lower limit
Pm15=1-interp1(grd,cdf,-15,'spline',0); %
return
