%%Plots weighted KDF stoichiometry distributions with fits previously acquired from cftool using histogram.

%% first load output.mat

%select appropriate channel
plottype=1; %(0) histogram; (1) kernel density estimate; (2) violin plot
fitflag = 0;
ch=1; %2;  %INPUT REQUIRED
showlinked=0;
showunlinked=0;
showtotal=1;
showPairedSingle=0;
maxstoich=1000;
maxnumtracks=10000;

%GW = params.stoichKDW;
GW = 0.6; %Gaussian width custom, default = 0.6-0.7  %don't change unless specifically looking for finer stoichiometry peaks at high exposure

%% Load output file

disp('Please load analysis output.mat');
[outputName,outputDir]=uigetfile('*.mat');
cd(outputDir)
load(outputName);

%% Get mean stats
try
    Ntracks = output.NlinkTrack+output.NunlinkTrack;
    StoichsMean = output.unlinkStoichsMean.*output.NunlinkTrack./Ntracks+output.linkStoichsMean.*output.NlinkTrack./Ntracks;
catch
    try %only unlinked
        Ntracks = output.NunlinkTrack;
        StoichsMean = output.unlinkStoichMean.*output.NunlinkTrack./Ntracks;
    catch
            error('No tracks found.')
    end
end


%% Plot S figures  (first load .sfit for linked and unlinked data)

%extract params from output.mat
if ~exist('params')
    if iscell(paramsA)
        params=paramsA{1,1};
    else
        params=paramsA;
    end
end

if ch==1
    dispcolorL = 'red';%'magenta';
    dispcolorU = 'blue';%'red';
    dispcolorT = 'black';
elseif ch==2
    dispcolorL = 'red';
    dispcolorU = 'blue';%'blue';
    dispcolorT = 'cyan';
end

%get histogram from which fits derived. Width of bins is same as kernel density width.
[stoichsU,edgesU] = histcounts(output.UnlinkedStoichsList{1,1},ceil(max(output.UnlinkedStoichsList{1,1})/GW));
binsU = edgesU(2:end)/2+edgesU(1:end-1)/2;

if showlinked && output.NlinkTrack(1) > 0
    [stoichsL,edgesL] = histcounts(output.LinkedStoichsList{1,1},ceil(max(output.LinkedStoichsList{1,1})/GW));
    binsL = edgesL(2:end)/2+edgesL(1:end-1)/2;
end

totalstoichs=[output.UnlinkedStoichsList{1,ch};output.LinkedStoichsList{1,ch}];
[stoichsT,edgesT] = histcounts(totalstoichs,ceil(max(totalstoichs)/GW));
binsT = edgesT(2:end)/2+edgesT(1:end-1)/2;

%% get violin plot of stoichiometry

if plottype == 2
    if showunlinked ==1
        if output.NlinkTrack == 0
            disp("No linked tracks found in this dataset.")
        end
        unlinkedstoichs=output.UnlinkedStoichsList{1,ch};
        if length(unlinkedstoichs)>maxnumtracks
            unlinkedstoichs=unlinkedstoichs(randperm(length(unlinkedstoichs)));
            unlinkedstoichs=unlinkedstoichs(1:maxnumtracks);
        end
        unlinkedstoichs=unlinkedstoichs(unlinkedstoichs<maxstoich);
        violinPlot(unlinkedstoichs(randperm(length(unlinkedstoichs))),ch);
    end
    if showlinked == 1
        if output.NlinkTrack == 0
            disp("No linked tracks found in this dataset.")
        else
            linkedstoichs=output.LinkedStoichsList{1,ch};
            if length(linkedstoichs)>maxnumtracks
                linkedstoichs=linkedstoichs(randperm(length(linkedstoichs)));
                linkedstoichs=linkedstoichs(1:maxnumtracks);
            end
            linkedstoichs=linkedstoichs(linkedstoichs<maxstoich);
            violinPlot(linkedstoichs(randperm(length(linkedstoichs))),ch);
        end
    end
    if showtotal == 1
        totalstoichs=totalstoichs(totalstoichs<maxstoich);
        violinPlot(totalstoichs(randperm(length(totalstoichs))),ch);
    end
    if showPairedSingle ==1
        pairedstoichs=output.SingleRbPairedStoichs(:,ch);
        pairedstoichs=pairedstoichs(pairedstoichs<maxstoich);
        pairedstoichs=pairedstoichs(randperm(length(pairedstoichs)));
        if length(pairedstoichs)>maxnumtracks
            pairedstoichs=pairedstoichs(1:maxnumtracks);
        end
        violinPlot(pairedstoichs,ch);
    end

    %pbaspect([1 1 1]);
    ylabel('Stoichiometry (molecules)')
    ylim([0,maxstoich])
    h=get(gca,'Children');
    set(h(4),'MarkerSize',1);
    box on
    
else
    %% get kernel density plot (not truncated at zero)
    if showunlinked ==1
        if output.NlinkTrack == 0
            disp("No linked tracks found in this dataset.")
        end
        if plottype ==1
            [plotU,KDFstoichsU,KDFxU]=KDFplotH(output.UnlinkedStoichsList{1,ch},GW);
            set(plotU,'Color',dispcolorU); hold on
        elseif plottype ==0
            htotalU = histogram(output.UnlinkedStoichsList{1,ch},'FaceColor',dispcolorU,'BinMethod','integers');
        end
    end
    
    if showlinked
        if output.NlinkTrack > 0
            if plottype ==1
                [plotL,KDFstoichsL,KDFxL]=KDFplotH(output.LinkedStoichsList{1,ch},GW);
                set(plotL,'Color',dispcolorL); hold on
            elseif plottype == 0
                htotalL = histogram(output.LinkedStoichsList{1,ch},'FaceColor',dispcolorL,'BinMethod','integers');
                hold on
            end
        else
            disp("No linked tracks found in this dataset.")
        end
    end
    
    if showtotal ==1
        if plottype ==1
            [plotU,KDFstoichsT,KDFxT]=KDFplotH(totalstoichs,GW);
            set(plotU,'Color',dispcolorT); hold on
        elseif plottype == 0
            htotalS = histogram(totalstoichs,'FaceColor',dispcolorT,'BinMethod','integers');
            hold on
        end
    end
    
    pbaspect([1 1 1]);
    ylabel('Track Probability Density')
    xlabel('Stoichiometry (molecules)')
    box on
    rangemax = max([binsU(end),binsU(end)])+10;  %set the axis limit to at least 6
    upsampledx=linspace(0,rangemax,500);  %upsample the number of points for plotting the fit components
    xlim([0,rangemax])
    
    %% FITTING
    if fitflag == 1
        
        cftool    %generate fits in cftool and save them as FitU.sfit, FitL.sfit
        
        
        %Unlikely to be fits specific for total data
        
        %Syntax is for Gaussian amplitudes in sequence first, then a single standard deviation variable, then the stoichiometry centres in sequence.
        
        % %LOAD FIT FOR LINKED DATA
        % if plottype==1
        % if ch==1
        %     FitL = Fit_stoichsDnaQlinkedKDF;             %INPUT REQUIRED
        %     FitU = Fit_stoichsDnaQunlinkedKDF;
        % elseif ch==2
        %     FitL = Fit_stoichsDnaTlinkedKDF;             %INPUT REQUIRED
        %     FitU = Fit_stoichsDnaTunlinkedKDF;
        % end
        % else
        %     if ch==1
        %     FitL = Fit_stoichsDnaQlinked;             %INPUT REQUIRED
        %     FitU = Fit_stoichsDnaQunlinked;
        % elseif ch==2
        %     FitL = Fit_stoichsDnaTlinked;             %INPUT REQUIRED
        %     FitU = Fit_stoichsDnaTunlinked;
        % end
        % end
        
        %unlinked
        if showunlinked == 1
            if ~exist('FitU')
                FitU=uigetfile('*.sfit');
            end
            fittedmodelU = FitU.AllFitdevsAndConfigs{1, 1}.Fitdev.Fit;
            CU = coeffvalues(fittedmodelU);  %extract fitted coefficients
            numcompsU = floor(length(CU)/2);
            if length(CU)/2 == numcompsU  %test if S width parameter doesn't exist
                S_existsU = 0;
                GWU = GW;
            else
                S_existsU = 1;
                GWU = CU(numcompsU+1);
            end
            for i = 1:numcompsU
                curveU(:,i) = CU(i)*exp(-(upsampledx-CU(i+S_existsU+numcompsU)).^2/(2*GWU^2));
            end
            totalU = sum(curveU(:,:),2);
            scalefactorU = max(totalU)/max(KDFstoichsU);
            
            for j = 1:numcompsU
                hold on; line(upsampledx,curveU(:,j)/scalefactorU,'Color','black','LineStyle','--');
            end
            hold on; line(upsampledx,totalU/scalefactorU,'Color',dispcolorU,'LineStyle','-');
            %hold on; bar(binsU,stoichsU/sum(stoichsU)/(binsU(2)-binsU(1)),'FaceColor','None','BarWidth',0.1);
            
        end
        
        %linked
        
        if showlinked == 1 && output.NlinkTrack > 0
            if ~exist('FitL')
                FitL=uigetfile('*.sfit');
            end
            fittedmodelL = FitL.AllFitdevsAndConfigs{1, 1}.Fitdev.Fit;
            CL = coeffvalues(fittedmodelL);  %extract fitted coefficients and define component distributions
            numcompsL = floor(length(CL)/2);
            if length(CL)/2 == numcompsL  %test if S width parameter doesn't exist
                S_existsL = 0;
                GWL = GW;
            else
                S_existsL = 1;
                GWL = CL(numcompsL+1);
            end
            for i = 1:numcompsL
                curveL(:,i) = CL(i)*exp(-(upsampledx-CL(i+S_existsL+numcompsL)).^2/(2*GWL^2));
            end
            totalL = sum(curveL(:,:),2);
            scalefactorL = max(totalL)/max(KDFstoichsL);
            
            for i = 1:numcompsL
                hold on; line(upsampledx,curveL(:,i)/scalefactorL,'Color','black','LineStyle','--');
            end
            hold on; line(upsampledx,totalL/scalefactorL,'Color',dispcolorL,'LineStyle','-');
            %hold on; bar(binsL,stoichsL/sum(stoichsL)/(binsL(2)-binsL(1)),'FaceColor','None','BarWidth',0.1);
        end
    end
end
