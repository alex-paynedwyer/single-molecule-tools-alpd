%% HOTPLOT
% creates a heatmap of scatterplot data
% X,Y 1D vectors containing data normally plotted with a scatter plot
% binSize sets size bin in a 2D histogram of the data
% Resolution sets the number of points to interpolate over to make a nice smooth map
% SetbinCentres use this to set the location of the bins, usually so you
% can plot on same axes (can use binCentres output) OPTIONAL
% heatmapdata is the matrix of intensity values
% binCentres marks the location of the bins


function [heatmapdata,PlotRange]=gausshotplot(X,Y,Resolution,PSFwidth,PlotRange)
if exist('PlotRange')==0
    PlotRange(1,1)=min(X);
    PlotRange(1,2)=max(X);
        PlotRange(2,1)=min(Y);
    PlotRange(2,2)=max(Y);
end

Intensity=1000;
image=zeros(Resolution,Resolution);
 [Xpos,Ypos] = meshgrid(1:Resolution,1:Resolution);
for i=1:length(X)
    Xpix=round(((X(i)-PlotRange(1,1))./(PlotRange(1,2)-PlotRange(1,1))).*Resolution);
    Ypix=round(((Y(i)-PlotRange(2,1))./(PlotRange(2,2)-PlotRange(2,1))).*Resolution);
    image=image+(Intensity/(2.*pi.*PSFwidth.^2))*exp(-(((Xpos-Xpix).^2)./(2.*PSFwidth^2)+((Ypos-Ypix).^2)./(2.*PSFwidth^2)));
end



heatmapdata=image;
surf(image','EdgeColor','none')
zlim([0.5,10000])
% xlim([min(min(Xq)),max(max(Xq))])
% ylim([min(min(Yq)),max(max(Yq))])set(gca,'XTickMode','manual');

view(2)
set(gca,'XTickMode','manual');
set(gca,'XTick',0:Resolution/5:Resolution)
set(gca,'XtickLabel',round(PlotRange(2,1):(PlotRange(2,2)-PlotRange(2,1))/5:PlotRange(2,2)));

set(gca,'YTickMode','manual');
set(gca,'YTick',0:Resolution/5:Resolution)
set(gca,'YtickLabel',round(PlotRange(1,1):(PlotRange(1,2)-PlotRange(1,1))/5:PlotRange(1,2)));

end