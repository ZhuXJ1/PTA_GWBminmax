function [AyrV,pdf,AyrM,AyrL,AyrU]=CompAmaxfun(m0,mf0,mfe,Nsim)
% Input: BH mass in log10 scale (m0), BHMF (mf0)
% uncertainty in BHMF (mfe), number of simulations (Nsim)
% Output: coordinates in Ayr^max (AyrV), the pdf of Ayr^max 
% and median and upper bound (AyrL, AyrU)

Aref=1.7228e-21; % constant in Eq. 10 of the paper
Amat=zeros(Nsim,1);
sjs=randn(1,Nsim);
%
zcut=2; % cutoff redshift
phidr = @(z) (1+z).^(-1/3);
zm13 = (1.1^(5/3))*integral(phidr,0,zcut);
m=10.^(m0);
for k=1:Nsim
    mfV=mf0+sjs(k).*mfe;
    % integrating over {M^(5/3)dn/dM}dM
    js1=trapz(m0,(10.^(mfV)).*(m.^(5/3)));
    Amat(k)=sqrt(zm13)*Aref.*sqrt(js1);
end
[pdf,grid]=kde1d(log10(Amat));
AyrV=10.^grid;
%
cdf=cumtrapz(grid,pdf);
logA95=interp1(cdf,grid,0.5);
AyrM=10^(logA95); % median Amax
logA95=interp1(cdf,grid,0.16);
AyrL=10^(logA95); % 1-sigma lower Amax
logA95=interp1(cdf,grid,0.84);
AyrU=10^(logA95); % 1-sigma upper Amax
return