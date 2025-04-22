%%Plots KDF Diffusivity distributions with Gamma fits previously acquired from cftool using histogram.

%% first load output.mat

%select appropriate channel
KDFflag=1; % 0 for histogram, 1 for kernel density plot
fitflag = 0; %set to 0 before fitting and 1 to plot an existing sfit that you loaded from cftool and renamed to FitU, FitL, FitT etc.
reflect=0;
ch = 1; %2;                         %INPUT REQUIRED
showlinked =0;  %choose one of the following
showunlinked =0;
showtotal =1;


%% Get mean stats

Ntracks = output.NlinkTrack+output.NunlinkTrack
NlinkTrack=output.NlinkTrack;
NunlinkTrack=output.NunlinkTrack;
try
diffKDW = params.diffKDW;
catch
diffKDW = paramsA{1,1}.diffKDW;
end
diffKDW=0.2;
try
    DiffsMean = output.unlinkDiffsMean.*output.NunlinkTrack./Ntracks+output.linkDiffsMean.*output.NlinkTrack./Ntracks
catch
    DiffsMean = output.unlinkDiffMean.*output.NunlinkTrack./Ntracks+output.linkDiffMean.*output.NlinkTrack./Ntracks
end

%% Plot D figures  (first load .sfit for linked and unlinked data)

%extract params from output.mat

if ~exist('params')
    if iscell(paramsA)
        params=paramsA{1,1};
        if iscell(params)
            params=params{1,1};
        end
    else
        params=paramsA;
    end
end

try
%gamma_width = diffKDW/2; %KDE width
gamma_width = sqrt(2)*max(output.unlinkDiffsMean(2,ch),diffKDW/sqrt(output.unlinkStoichsMean(1,ch))); %KDE width linked to stoichiometry
%gamma_width = 1; % manual KDE width
hist_width = max(output.unlinkDiffsMean(2,ch),diffKDW/sqrt(output.unlinkStoichsMean(1,ch)));
%hist_width = 0.09;
diffmin = max(sqrt(2)*output.unlinkDiffsMean(2,ch),sqrt(2)*diffKDW/sqrt(output.unlinkStoichsMean(1,ch))); %cutoff for nonphysical D coefficient (e.g. negative gradient in getDiffusion3)
catch
%gamma_width = diffKDW/2; %KDE width
gamma_width = max(output.unlinkDiffsMean(ch,2),diffKDW/sqrt(output.unlinkStoichsMean(ch,1))); %KDE width linked to stoichiometry
%gamma_width = 0.08; % manual KDE width
%hist_width = max(output.unlinkDiffsMean(ch,2),diffKDW/sqrt(output.unlinkStoichsMean(ch,1)));
hist_width = 0.006;
diffmin = max(sqrt(2)*output.unlinkDiffsMean(ch,2),sqrt(2)*diffKDW/sqrt(output.unlinkStoichsMean(ch,1))); %cutoff for nonphysical D coefficient (e.g. negative gradient in getDiffusion3)
end    

diffmin=0;%.0001;

if ch==1
    dispcolorL = 'red';%'magenta';
    dispcolorU = 'blue';
    dispcolorT = 'magenta';
elseif ch==2
    dispcolorL = 'blue';%'green';
    dispcolorU = 'red';%'blue';
    dispcolorT = 'cyan';
end

%get histogram from which fits derived.
[diffsU,edgesU] = histcounts(output.UnlinkedDiffsList{1,ch}(output.UnlinkedDiffsList{1,ch}>diffmin),ceil(max(output.UnlinkedDiffsList{1,ch}(output.UnlinkedDiffsList{1,ch}>diffmin))/hist_width));
binsU = edgesU(2:end)/2+edgesU(1:end-1)/2;

if showlinked %&& output.NlinkTrack 
[diffsL,edgesL] = histcounts(output.LinkedDiffsList{1,ch}(output.LinkedDiffsList{1,ch}>diffmin),ceil(max(output.LinkedDiffsList{1,ch}(output.LinkedDiffsList{1,ch}>diffmin))/hist_width));
binsL = edgesL(2:end)/2+edgesL(1:end-1)/2;
end

totaldiffs=[output.UnlinkedDiffsList{1,ch}(output.UnlinkedDiffsList{1,ch}>diffmin);output.LinkedDiffsList{1,ch}(output.LinkedDiffsList{1,ch}>diffmin)];
[diffsT,edgesT] = histcounts(totaldiffs,ceil(max(totaldiffs)/hist_width));
binsT = edgesT(2:end)/2+edgesT(1:end-1)/2;

%get kernel density plot (not truncated at zero)

if showunlinked ==1
    if KDFflag ==1
        [plotU,KDFdiffsU,KDFxU]=KDFplotH(output.UnlinkedDiffsList{1,ch}(output.UnlinkedDiffsList{1,ch}>diffmin),gamma_width);
        set(plotU,'Color',dispcolorU); hold on
    else
        hunlinkedD = histogram(output.UnlinkedDiffsList{1,ch}(output.UnlinkedDiffsList{1,ch}>diffmin),'FaceColor',dispcolorU,'BinWidth',hist_width);
        hold on
    end
end

if showlinked %&& NlinkTrack > 0 
    if KDFflag ==1
[plotL,KDFdiffsL,KDFxL]=KDFplotH(output.LinkedDiffsList{1,ch}(output.LinkedDiffsList{1,ch}>diffmin),gamma_width);
   set(plotL,'Color',dispcolorL); hold on
    else
      hlinkedD = histogram(output.LinkedDiffsList{1,ch}(output.LinkedDiffsList{1,ch}>diffmin),'FaceColor',dispcolorL,'BinWidth',hist_width);
        hold on   
        
    end
end

if showtotal ==1
    if KDFflag ==1
        gamma_width=0.1;
    [plotU,KDFdiffsT,KDFxT]=KDFplotH(totaldiffs,gamma_width);
    set(plotU,'Color',dispcolorT); hold on
    else
        htotalD = histogram(totaldiffs,'FaceColor',dispcolorT,'BinWidth',hist_width);
        hold on
    end
end
pbaspect([1 1 1]);
xlabel('Diffusion Coefficient ({\mu}m^{2}/s)')
ylabel('Probability Density')
box on

rangemax = max([binsU(end),binsU(end),0.1]);  %set the axis limit to at least 0.1 sq micron/s
%rangemax = max([binsU(end),binsU(end),1]);  %set the axis limit to at least 1 sq micron/s
%rangemax = 1.4;
upsampledx=linspace(0,rangemax,500);  %upsample the number of points for plotting the fit components
xlim([0,rangemax])
%ylim([0,2])

if fitflag == 1

%cftool    %generate fits and save them
% Gamma fit is (copy to cftool custom equation): y = a*gampdf(x,d,e/(d-1))+b*gampdf(x,d,f/(d-1))+c*gampdf(x,d,g/(d-1))
% Syntax is for Gamma amplitudes in sequence, then a common shape variable, then the modal diffusities in sequence.

% %LOAD FIT FOR LINKED DATA
    
% if KDFflag==1
% if ch==1
%     FitL = Fit_diffsDnaQlinkedKDF;             %INPUT REQUIRED
%     FitU = Fit_diffsDnaQunlinkedKDF;
% elseif ch==2
%     FitL = Fit_diffsDnaTlinkedKDF;             %INPUT REQUIRED
%     FitU = Fit_diffsDnaTunlinkedKDF;
% end
% else
%     if ch==1
%     FitL = Fit_diffsDnaQlinked;             %INPUT REQUIRED
%     FitU = Fit_diffsDnaQunlinked;
% elseif ch==2
%     FitL = Fit_diffsDnaTlinked;             %INPUT REQUIRED
%     FitU = Fit_diffsDnaTunlinked;
% end
% end

%unlinked
if showunlinked == 1

fittedmodelU = FitU.AllFitdevsAndConfigs{1, 1}.Fitdev.Fit;
CU = coeffvalues(fittedmodelU);  %extract fitted coefficients
curve1U = CU(1)*gampdf(upsampledx,CU(4),CU(5)/(CU(4)-1)); %define component distributions
curve2U = CU(2)*gampdf(upsampledx,CU(4),CU(6)/(CU(4)-1));
curve3U = CU(3)*gampdf(upsampledx,CU(4),CU(7)/(CU(4)-1));
totalU = curve1U + curve2U + curve3U;
scalefactorU = CU(1)+CU(2)+CU(3);
%scalefactorU = max(totalU)/max(diffsU);
%scalefactorU = max(totalU)/max(KDFdiffsU);

hold on; bar(binsU,diffsU/sum(diffsU)/(binsU(2)-binsU(1)),'LineStyle','None','FaceColor',[0.6,0.6,0.6],'BarWidth',1);
hold on; line(upsampledx,curve1U/scalefactorU,'Color','blue','LineStyle','-','LineWidth',2);
hold on; line(upsampledx,curve2U/scalefactorU,'Color','red','LineStyle','-','LineWidth',2);
hold on; line(upsampledx,curve3U/scalefactorU,'Color','green','LineStyle','-','LineWidth',2);
hold on; line(upsampledx,(curve1U+curve2U+curve3U)/scalefactorU,'Color','black','LineStyle','--','LineWidth',2);

end

if showlinked & NlinkTrack > 0 

fittedmodelL = FitL.AllFitdevsAndConfigs{1, 1}.Fitdev.Fit;
CL = coeffvalues(fittedmodelL);  %extract fitted coefficients and define component distributions

curve1L = CL(1)*gampdf(upsampledx,CL(4),CL(5)/(CL(4)-1));
curve2L = CL(2)*gampdf(upsampledx,CL(4),CL(6)/(CL(4)-1));
curve3L = CL(3)*gampdf(upsampledx,CL(4),CL(7)/(CL(4)-1));
totalL = curve1L + curve2L + curve3L;
scalefactorL = CL(1)+CL(2)+CL(3);
%scalefactorL = max(totalL)/max(diffsL);
%scalefactorL = max(totalL)/max(KDFdiffsL);

hold on; bar(binsL,diffsL/sum(diffsL)/(binsL(2)-binsL(1)),'LineStyle','None','FaceColor',[0.6,0.6,0.6],'BarWidth',1);
hold on; line(upsampledx,curve1L/scalefactorL,'Color','blue','LineStyle','-','LineWidth',2);
hold on; line(upsampledx,curve2L/scalefactorL,'Color','red','LineStyle','-','LineWidth',2);
hold on; line(upsampledx,curve3L/scalefactorL,'Color','green','LineStyle','-','LineWidth',2);
hold on; line(upsampledx,(curve1L+curve2L+curve3L)/scalefactorL,'Color','black','LineStyle','--','LineWidth',2);

elseif showlinked == 1 & NlinkTrack < 0
disp('no linked tracks found')
end

if showtotal == 1

fittedmodelT = FitT.AllFitdevsAndConfigs{1, 1}.Fitdev.Fit;
CT = coeffvalues(fittedmodelT);  %extract fitted coefficients
curve1T = CT(1)*gampdf(upsampledx,CT(4),CT(5)/(CT(4)-1)); %define component distributions
curve2T = CT(2)*gampdf(upsampledx,CT(4),CT(6)/(CT(4)-1));
curve3T = CT(3)*gampdf(upsampledx,CT(4),CT(7)/(CT(4)-1));
totalT = curve1T + curve2T + curve3T;
scalefactorT = CT(1)+CT(2)+CT(3);
%scalefactorT = max(totalT)/max(diffsT);
%scalefactorU = max(totalT)/max(KDFdiffsT);

hold on; bar(binsT,diffsT/sum(diffsT)/(binsT(2)-binsT(1)),'LineStyle','None','FaceColor',[0.6,0.6,0.6],'BarWidth',1);
hold on; line(upsampledx,curve1T/scalefactorT,'Color','blue','LineStyle','-','LineWidth',2);
hold on; line(upsampledx,curve2T/scalefactorT,'Color','red','LineStyle','-','LineWidth',2);
hold on; line(upsampledx,curve3T/scalefactorT,'Color','green','LineStyle','-','LineWidth',2);
hold on; line(upsampledx,(curve1T+curve2T+curve3T)/scalefactorT,'Color','black','LineStyle','--','LineWidth',2);

end

end