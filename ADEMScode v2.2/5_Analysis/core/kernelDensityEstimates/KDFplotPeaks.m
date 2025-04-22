function [Dens, x,h]=KDFplotPeaks(data,bandwidthval)
color=rand(1,3);
if exist('bandwidthval')==0
    [Dens,x] = ksdensity(data,'npoints',10000);
    h=plot(x,Dens,'color',color);
else
    [Dens,x] = ksdensity(data,'npoints',10000,'bandwidth',bandwidthval);
    h=plot(x,Dens,'color',color);
   
end
 [pks, locs]=findpeaks(Dens);
 hold on
 scatter(x(locs),pks)
 %text(x(locs),pks,num2str(round(x(locs))))
 for t=1:length(locs)
     p=locs(t);
     text(x(p),pks(t)+max(pks)*0.1,num2str(x(p),3))
 end
end