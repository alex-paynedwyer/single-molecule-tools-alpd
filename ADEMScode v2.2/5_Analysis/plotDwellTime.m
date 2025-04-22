%%Plots dwell times of linked tracks with exponential decay fits previously acquired from cftool using histogram.

%load output.mat
%load .sfit for linked and unlinked data

%extract params from output.mat
if ~exist('params')
    if iscell(paramsB)
        params=paramsB{1,1};
    else
        params=paramsB;
    end
end

%select appropriate channel
ch = 1; %2;                         %INPUT REQUIRED
frametime = params.frameTime;
binwidth = round(1000*frametime,1);

if ch==1
    spotset = output.AllSpots1linked;
    trackset = output.trackArrayCh1;
    
    dispcolor = 'black';%'magenta';

elseif ch==2
    spotset = output.AllSpots2linked;
    trackset = output.trackArrayCh2;
       
    dispcolor = 'black';%'green';
end

%get track lengths if recorded
%lentracksU = trackset(trackset(:,10)==0,14);
%lentracksL = trackset(trackset(:,10)>0,14);
%lentracksT = trackset(:,14);


%get histogram from which fits derived.
%dwellset = trackset(trackset(:,10)>0,9)*params.frameTime*1000;
dwellset = trackset(:,9)*frametime*1000;
[dwellcounts,dwelledges] = histcounts(dwellset,ceil(max(dwellset)),'BinWidth',binwidth);
dwellbins = dwelledges(2:end)/2+dwelledges(1:end-1)/2;

disp('Mean/SEM colocalisation times')
[mean(dwellset(dwellset>0));std(dwellset(dwellset>0))/sqrt(length(dwellset(dwellset>0)))]

%cftool    %generate zero-intercept exponential fits f=a*exp(-b*x) and save them
%Syntax is for a,b (amplitude then exponent).

%plot together 

linked_dwell = histogram(dwellset,'FaceColor',dispcolor,'BinWidth',binwidth);
pbaspect([1 1 1]);
xlabel('Dwell time (ms)')
ylabel('Number of detected foci')
box on
hold on

rangemax = 80;%max([dwellbins(end),4000*frametime]);  %set the axis limit to at least 4 frames
upsampledx=linspace(0,rangemax,500);  %upsample the number of points for plotting the fit components
xlim([0,rangemax+binwidth/2])
%ylim([0,])

%cftool    %generate fits and save them

% %LOAD FIT FOR LINKED DWELL TIMES
%fittedmodel = fitresult;
fittedmodel = dwelltime.AllFitdevsAndConfigs{1, 1}.Fitdev.Fit;
CU = coeffvalues(fittedmodel);  %extract fitted coefficients
curve1U = CU(1)*exp(-upsampledx/CU(2));%define component distributions
hold on; line(upsampledx,curve1U,'Color','blue','LineStyle','-','LineWidth',2);
