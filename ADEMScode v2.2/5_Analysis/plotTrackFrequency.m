%PlotTrackFrequency
figure
if ~exist('params')
    if iscell(paramsA)
        params=paramsA{1,1};
    else
        params=paramsA;
    end
end

%select appropriate channel
ch =1; %1;                         %INPUT REQUIRED

showlinked =0;
showunlinked =0;
showtotal =1;

if ch==1
    dispcolorL = 'red';
    dispcolorU = 'blue';
    dispcolorT = 'magenta';
elseif ch==2
    dispcolorL = 'red';
    dispcolorU = 'blue';
    dispcolorT = 'green';
end

NCells=output.NCells(ch);

[trajsU,edgesU] = histcounts(output.UnlinkedTrajList{1,ch},NCells,'BinMethod','integers');
binsU = edgesU(2:end)/2+edgesU(1:end-1)/2;

[trajsL,edgesL] = histcounts(output.LinkedTrajList{1,ch},NCells,'BinMethod','integers');
binsL = edgesL(2:end)/2+edgesL(1:end-1)/2;

totaltrajs=[output.UnlinkedTrajList{1,ch};output.LinkedTrajList{1,ch}];
[trajsT,edgesT] = histcounts(totaltrajs,NCells,'BinMethod','integers');
binsT = edgesT(2:end)/2+edgesT(1:end-1)/2;

if showunlinked ==1
    hunlinked = histogram(trajsU,'FaceColor',dispcolorU);
    hold on
end
if showlinked ==1
    hlinked = histogram(trajsL,'FaceColor',dispcolorL);
    hold on
end
if showtotal ==1
    htotal = histogram(trajsT,'FaceColor',dispcolorT);
    hold on
end

xlabel('Tracks per cell')
ylabel('Frequency (cells)')
axis tight
pbaspect([1 1 1]);
