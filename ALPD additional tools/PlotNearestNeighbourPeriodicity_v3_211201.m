%Plot periodicity based on nearest neighbour peak distances

%load output file first

%% Calculate nearest neighbour spacing
figure; hold on
maxspacing = 60; %maximum credible interval
basePeakGW = 0.6; %molecular precision to display/smooth intervals
ch = 2;
colour='g'; %'g' or 'r'
total = 0;
linked = 1;

GW =0.6;%params.stoichKDW; %molecular precision to find intervals
rangemax = 60;
Isingle_interp_timepts = 3;  %number of points over which stoichiometry is back-interpolated for a single molecule
interp_timepts = 4;%params.uppertracklimit;  %mean number of points over which stoichiometry is back-interpolated for the current data
plotstoichs=0;
plotviolin=0;
oldcode=0;
selfweighted=0;
simulate=0;
correctionfactor=1;%1 for all labelled FPs (no endogenous protein)

%% back compatibility with old code that doesn't have lists in output file
if simulate
else
    if ch==1
    try
        Isingle = params.IsingleCh1;
    catch
        if oldcode==1
        if colour == 'r'
            Isingle = params.IsingleR;
        elseif colour == 'g'
            Isingle = params.IsingleG;
        end
        end
    end
elseif ch==2
    try
        Isingle = params.IsingleCh2;
    catch
        if oldcode==1
        Isingle = params.IsingleG;
        end
    end
end
end

if simulate==1
    periodgt=2;
    mu=21;
    listlen = 460;
    noiseamp=GW;
    [stoichgt,stoichs]=SimulatePoissonStoichiometry(periodgt,mu,GW,listlen,1);figure;
else
    if total==0
        if linked==1
            if oldcode==1
                if ch==1
                    stoichs=output.trackArrayCh1(output.trackArrayCh1(:,8)>0,2)/Isingle;
                elseif ch==2
                    stoichs=output.trackArrayCh1(output.trackArrayCh2(:,8)>0,2)/Isingle;
                end
            else
                stoichs=output.LinkedStoichsList{1,ch};
            end
        elseif linked==0
            if oldcode==1
                if ch==1
                    stoichs=output.trackArrayCh1(output.trackArrayCh1(:,8)==0,2)/Isingle;
                elseif ch==2
                    stoichs=output.trackArrayCh1(output.trackArrayCh2(:,8)==0,2)/Isingle;
                end
            else
                stoichs=output.UnlinkedStoichsList{1,ch};
            end
        end
    elseif total==1
        if oldcode ==1
            if ch==1
                stoichs=output.trackArrayCh1(:,2)/Isingle;
            elseif ch==2
                stoichs=output.trackArrayCh2(:,2)/Isingle;
            end
        else
            stoichs=[output.LinkedStoichsList{1,ch};output.UnlinkedStoichsList{1,ch}];
        end        
    end
end

%%
stoichs=stoichs*correctionfactor;

if selfweighted
    [plotP,KDFpersL,KDFxP]=KDFplotH(stoichs,GW,0,1);
else
    [plotP,KDFpersL,KDFxP]=KDFplotH(stoichs,GW);
end
if plotstoichs == 0
    close gcf
end
[pks,locs] = findpeaks(KDFpersL);
peaklocs=KDFxP(locs);
peaklocs=[0,peaklocs];
spacing=(peaklocs(2:end)-peaklocs(1:end-1));
spacing=spacing(spacing<maxspacing); %remove spaces higher than this as spurious multiples;
PeakGW = basePeakGW*sqrt(2*mean(stoichs)/length(locs)/(interp_timepts/Isingle_interp_timepts))/sqrt(2);
figure
if selfweighted==1
    [plotS,KDFperS,KDFxS]=KDFplotH(spacing(spacing>GW),PeakGW,0,1);
else
    [plotS,KDFperS,KDFxS]=KDFplotH(spacing(spacing>GW),PeakGW,0,0);
end

[pksmode,locsmode]=findpeaks(KDFperS);
modalloc=locsmode(pksmode==max(pksmode));
modalcentre=KDFxS(modalloc(1))
nearmodalspacing =spacing(abs(spacing-modalcentre)<(GW*2));
mean_mode=mean(nearmodalspacing);
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