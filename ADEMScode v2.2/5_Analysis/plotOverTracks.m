function plotOverTracks(selectinput,spotsBaseline,MinL,MaxS,MinS,MinStep,Isingle)
if nargin<7
    selectinput=1;%#ok<*NOPRT> %heuristics for filtering the overtracked traces from a BASELINE.mat file
    
    Isingle=126; % nominal characteristic molecular brightness in ADU counts (don't make this too small to start with)
    MinL=3; % only keep tracks with this number of frames or greater
    AbsoMaxS=200;  % highest max intensity (as multiple of nominal Isingle)
    AbsoMinS=-0.9;  % lowest min intensity
    MaxS=0.7;  % lowest max intensity
    MinS=50;  % highest min intensity
    MinStep=1;  % minimum number of Isingle steps (careful to keep smaller than MaxS-MinS); 
  
    NstartCalc=1;   %initial track that is considered. Use with the MaxTrCalc to batch process tracks and whittle down to good examples
    MaxTrCalc=10000; %for time/memory's sake, limit the number of tracks that can be shown
end

if selectinput==1
    [DataFile,DataDir]=uigetfile('*.mat','Select the correct overtracking data preferably in BASELINE.mat');
    cd(DataDir);
    load(DataFile);
end

%% Plot the CK-filtered tracks which show single molecule steps

XXall=[];
spotindex=0;

baselinespots=spotsBaselineCh1;
%Isingle=params.IsingleCh1;
%frameTime=params.frameTime;
frameTime=1;
plotTraces=1;

tracknos=unique(baselinespots(:,10));
tracknos=tracknos(tracknos>0);kappalist=[];
Ntracks4calc=length(tracknos);
ncalc=min(Ntracks4calc-NstartCalc,MaxTrCalc-NstartCalc);

for alpha =NstartCalc:NstartCalc+ncalc
    kappa=tracknos(alpha);
    spots=sortrows(baselinespots(baselinespots(:,10)==kappa,:),9);
    %select only rows showing stepwise photobleaching of single molecule steps
    if min(spots(:,5))<MinS*Isingle && max(spots(:,5))>MaxS*Isingle && max(spots(:,5))<AbsoMaxS*Isingle && min(spots(:,5))>AbsoMinS*Isingle && length(spots(:,9))>MinL && isscalar(max(spots(:,9))) && isscalar(min(spots(:,9))) && (spots(spots(:,9)==min(spots(:,9)),5))>(spots(spots(:,9)==max(spots(:,9)),5)+Isingle*MinStep)
        spotindex=spotindex+1;
        [XX,~,~,~,~,~]=ckfiltb2original(spots(:,5),2,50);
        kappalist=[kappalist,kappa];
        N=length(spots(:,5));
        if plotTraces==1
			hold on;
			s=scatter((1:N)*frameTime,spots(:,5));
			s.Marker='o';
			s.SizeData=6;
			s.MarkerFaceColor='flat';
			line((1:N)*frameTime,XX,'Color',[0.5,0.5,0.5,0.5]);
		end
        %XXall=[XXall;XX];
        XXall=[XXall;XX-mean(XX(end-3:end))]; %realign to zero background
        %hold on;
    end
end

xlabel('Cumulative exposure (ms)');
ylabel('Track intensity (counts)');
end
