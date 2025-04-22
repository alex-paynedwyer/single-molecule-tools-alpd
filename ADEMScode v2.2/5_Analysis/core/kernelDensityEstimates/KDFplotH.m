% KDFplot which outputs a figure handle

function [h,Dens, x]=KDFplotH(data,bandwidthval,reflect,selfweight)
colorVal=rand(1,3);
if nargin<4
    selfweight=0;
end
if nargin<3
    if nargin==2
        [Dens,x] = ksdensity(data,'npoints',10000,'bandwidth',bandwidthval);
    elseif nargin==1
        [Dens,x] = ksdensity(data,'npoints',10000);
    end
else
    if reflect==0
        [Dens,x] = ksdensity(data,'npoints',10000,'bandwidth',bandwidthval);
    else
        [Dens,x] = ksdensity(data,'npoints',10000,'bandwidth',bandwidthval,'support',[-0.0000000001,2*max(data)],'BoundaryCorrection','reflection');
    end
end
if selfweight==1
h=plot(x,(Dens.*x)./(sum(Dens.*x))/mean(x(2:end)-x(1:end-1)),'color',colorVal,'LineWidth',2);
else
h=plot(x,Dens,'color',colorVal,'LineWidth',2);
end
end