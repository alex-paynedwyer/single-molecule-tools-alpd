%Plot periodicity based on nearest neighbour peak distances

%load output file first

%% Calculate nearest neighbour spacing
figure; hold on
maxspacing = 10; %maximum credible interval
stoichmax = 50; %max stoichiometry
basePeakGW = 0.3; %molecular precision to display/smooth intervals
ch = 1;
colour='g'; %'g' or 'r'
total = 0;
linked = 0;

GW =0.6;%params.stoichKDW; %molecular precision to find intervals
rangemax = 10;
manualIsingle = 126; %manual input if cannot be found from params
Isingle_interp_timepts = 3;  %number of points over which stoichiometry is back-interpolated for a single molecule
interp_timepts = 3;%params.uppertracklimit;  %mean number of points over which stoichiometry is back-interpolated for the current data
plotstoichs=1;
plotviolin=0;
oldcode=0;
simulate=0;
correctionfactor=1;%1 for all labelled FPs (no endogenous protein)

%% back compatibility with old code that doesn't have lists in output file
if simulate
else
    if ch==1
        try
            Isingle = params.IsingleCh1;
        catch
            Isingle = paramsA.IsingleCh1;
        end
    elseif ch==2
        try
            Isingle = params.IsingleCh2;
        catch
            Isingle = paramsA.IsingleCh2;
        end
    end
end
if ~exist('Isingle')
    Isingle=manualIsingle;
end

if simulate==1
    periodgt=2;%20*rand()+5;
    mu=9.2;%9*rand();
    listlen = 15*mu;
    noiseamp=50*GW;
    [stoichgt,stoichs]=SimulatePoissonStoichiometry(periodgt,mu,GW,listlen,1);figure;
else
    if total==0
        if linked==1
            stoichs=output.LinkedStoichsList{1,ch};
        elseif linked==0
            stoichs=output.UnlinkedStoichsList{1,ch};
        end
    elseif total==1
        stoichs=[output.LinkedStoichsList{1,ch};output.UnlinkedStoichsList{1,ch}];
    end
end

%%
stoichs=stoichs(stoichs<stoichmax);
stoichs=stoichs*correctionfactor;

[plotP,KDFpersL,KDFxP]=KDFplotH(stoichs,GW);
if plotstoichs == 0
    close gcf
end
[pks,locs] = findpeaks(KDFpersL);
peaklocs=KDFxP(locs);
peaklocs=[0,peaklocs];
spacing=(peaklocs(2:end)-peaklocs(1:end-1));
minspacing=max(GW,0.5);
spacing=spacing(spacing<maxspacing); %remove spaces higher than this as spurious multiples;
pks=pks(spacing<maxspacing);
spacing=spacing(spacing>minspacing); %remove spaces smaller than GW as spurious noise;
pks=pks(spacing>minspacing);
PeakGW = basePeakGW*sqrt(mean(stoichs)/length(locs)/(interp_timepts/Isingle))/2;
figure
[KDFperS, KDFxS]=KDFplotWeight(spacing,pks,PeakGW);

[pksmode,locsmode]=findpeaks(KDFperS);
modalloc=locsmode(pksmode==max(pksmode));
modalcentre=KDFxS(modalloc(1));
nearmodalspacing =spacing(abs(spacing-modalcentre)<(GW*2));
mean_mode=mean(nearmodalspacing)
if simulate==1
    sem_mode=1.5*std(nearmodalspacing)/sqrt(length(nearmodalspacing)-1)
else
    sem_mode=1.5*std(nearmodalspacing)/sqrt(length(nearmodalspacing)-1)
end

pbaspect([1 1 1]);
ylabel('Peak fraction')
xlabel('Stoichiometry peak-to-peak interval (molecules)')
box on
%rangemax = max(spacing)+2*GW;
xlim([0,rangemax])

if plotviolin==1
    violinPlot(spacing',ch);
end