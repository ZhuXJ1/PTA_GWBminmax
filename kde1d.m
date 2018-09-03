function [pdf,grid]=kde1d(X)
% simple kernel density estimation code for 1-d data
% using Gaussian kernel and a rule-of-thumb kernel width

n=length(X);
MAX=max(X,[],1);MIN=min(X,[],1);scaling=MAX-MIN;
MAX=MAX+scaling/10;MIN=MIN-scaling/10;scaling=MAX-MIN;
grid=(MIN:scaling/(2^8-1):MAX)';

sigma=std(X);
hwidth=1.06*sigma*(n^(-1/5)); % rule of thumb
pdf=zeros(length(grid),1);
for j = 1:length(grid)
    gau1 = 1/(sqrt(2*pi)*hwidth).*exp(-0.5.*((grid(j)-X).^2)./(hwidth.^2)); % hist
    pdf(j) = sum(gau1)/n;
end
jnm=trapz(grid,pdf);
pdf=pdf./jnm;
end