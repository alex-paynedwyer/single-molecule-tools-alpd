%Simulate multiple rounds of periodicity calculation based on nearest neighbour peak distances
function [mean_mode,sem_mode,KDFxS,KDFperS,nintervals]=SimulatePoissonPeriodicity(periodgt,mu,listlen,GW,maxstoichs,maxspacing,basePeakGW,corrfactor,plotflag)
% periodgt=3; %ground truth periodicity
% mu=50;  %mean stoichiometry
% listlen = 125;  %number of tracks in each iteration
% 
% GW = 0.6; %molecular precision to find intervals in stoichiometry
% maxspacing = 30; %maximum credible interval in stoichiometry
% basePeakGW = 0.3; %molecular precision to display/smooth the intervals (upscaled by the root mean stoichiometry)
% corrfactor=1; %stoichiometry scaling factor to mimic partial labelling

noiseamp=GW;  %the noise on each track's stoichiometry is at least the molecular precision
Isingle_interp_timepts = 3;  %number of points over which stoichiometry is back-interpolated for a single molecule
interp_timepts = 4;%params.uppertracklimit;  %mean number of points over which stoichiometry is back-interpolated for the current data

[stoichgt,stoichs]=SimulatePoissonStoichiometry(periodgt,mu,GW,listlen,plotflag);
stoichs=stoichs*corrfactor;
stoichs=stoichs(stoichs<maxstoichs);
[KDFpersL,KDFxP] = ksdensity(stoichs,'npoints',10000,'bandwidth',GW);

[pks,locs] = findpeaks(KDFpersL);
peaklocs=KDFxP(locs);
peaklocs=[0,peaklocs];
spacing=(peaklocs(2:end)-peaklocs(1:end-1));
minspacing=max(GW,0.5);
spacing=spacing(spacing<maxspacing); %remove spaces higher than this as spurious multiples;
pks=pks(spacing<maxspacing);
spacing=spacing(spacing>minspacing); %remove spaces smaller than GW as spurious noise;
pks=pks(spacing>minspacing);
nintervals=length(spacing);
PeakGW = basePeakGW*sqrt(mean(stoichs)/length(locs)/(interp_timepts/Isingle))/2;
[KDFperS, KDFxS]=KDFplotWeight(spacing,pks,PeakGW);

[pksmode,locsmode]=findpeaks(KDFperS);
modalloc=locsmode(pksmode==max(pksmode));
modalcentre=KDFxS(modalloc(1));
nearmodalspacing =spacing(abs(spacing-modalcentre)<GW*2);
mean_mode=mean(nearmodalspacing);
sem_mode=1.5*std(nearmodalspacing)/sqrt(length(nearmodalspacing)-1);

% pbaspect([1 1 1]);
% ylabel('Peak fraction')
% xlabel('Stoichiometry peak-to-peak interval (molecules)')
% box on
% rangemax = 10;
% %rangemax = max(spacing)+2*GW;
% xlim([0,rangemax])
% violinPlot(spacing',ch)
end