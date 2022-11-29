function v = violinPlot(xdatalist,colourchannel,xNames,yLabel,smooth)
if nargin <5        %defaults
    smooth = 1;   
    if nargin <4
        xNames = "";
        yLabel = "";
        if nargin <2
            colourchannel = 0;
        end
    end
end

%plot normal violin and box for strictly positive continuous data, e.g. copy, stoichiometry, number of tracks

figure
set(gcf,'Position',[1000 300 400 600])
distributionPlot(xdatalist,'divFactor',2,'histOpt',1,'addSpread',1,'showMM',4,'addBoxes',1,'yLabel',yLabel,'xNames',xNames);

if colourchannel==0
    patchcolor = [0.8 0.8 0.8];  % monochrome
    spotcolor = [0.4 0.4 0.4];
elseif colourchannel==1
    patchcolor = [1 0.8 1];      % magenta
    spotcolor = [0.8 0.5 0.8];
elseif colourchannel==2
    patchcolor = [0.7 0.9 0.7];    % green
    spotcolor = [0.3 0.7 0.3];
end

items = get(gca,'Children');
meanerrorbar = items(1);
set(meanerrorbar,'LineWidth',2.5,'Color','black');
meanmarker = items(2);
set(meanmarker,'LineWidth',2,'Color','black');
boxlist = items(3); boxlines = get(boxlist,'Children'); 

outliers = boxlines(1);
medianline = boxlines(2);
box = boxlines(3);
lowerquartile = boxlines(4);
upperquartile = boxlines(5);
lowerwhisker = boxlines(6);
upperwhisker = boxlines(7);

for i=1:length(boxlines)
set(boxlines(i),'Color','black','LineStyle','-','LineWidth',1);
end

delete(outliers)

spread = items(4); set(spread,'MarkerSize',10,'Color',spotcolor)
patch = items(5); set(patch,'FaceColor',patchcolor)

pbaspect([0.5 1 1]);
ylim([0,max(xdatalist)*1.1])