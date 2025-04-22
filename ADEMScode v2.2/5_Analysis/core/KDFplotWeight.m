function [Dens, x]=KDFplotWeight(data,weight,bandwidth)
% plots a weighted KDF based on multiplying by weight
if exist('bandwidth')==0
    [Dens,x] = ksdensity(data,'npoints',10000,'Weights',weight);
    h=plot(x,Dens,'LineWidth',2);
else
    [Dens,x] = ksdensity(data,'npoints',10000,'Weights',weight,'bandwidth',bandwidth);
    h=plot(x,Dens,'LineWidth',2);
end
end