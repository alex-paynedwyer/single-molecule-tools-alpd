function [Dens, x]=KDF(data,bandwidthval)
color=rand(1,3);
if exist('bandwidthval')==0
    [Dens,x] = ksdensity(data,'npoints',10000);
   % plot(x,Dens,'color',color)
else
    [Dens,x] = ksdensity(data,'npoints',10000,'bandwidth',bandwidthval);
   % plot(x,Dens,'color',color)
    
end
end