%% Runs a Chung-Kennedy filter over the data
% See ckfitb2original for full details
% Set bandwidth to roughly the on-time of a fluorophore
% weightval is between 1-100, start at 50
% minStep and maxStep filter the steps - choose Isingle to fall between these values.

%e.g. output=CKall_ALPD(spotsBaselineCh1,5,100,10000,50,1)

function filteredI=CKfilterBaseline(spots,bandwidth,minStep,maxStep,weightVal,showOutput)
AllXX=[];
i=0; mintraces=350; maxtraces=500;
for q=1:max(spots(:,10))
    lengthSpot=size(spots(spots(:,10)==q,9),1);
    SpotI=spots(spots(:,10)==q,5);
    %    SpotStoich(i)=mean(SpotI(1:end))/2500;

    try
        [XX,TX,DX,SD,DSD,XPRE]=ckfiltb2original(spots(spots(:,10)==q,5),bandwidth,weightVal);
        AllXX=[AllXX,XX'];
       
        %
        if showOutput==1 && (max(SpotI)-min(SpotI))>(minStep)  && (max(SpotI)-min(SpotI))<(maxStep)   %try to plot only single molecule steps
            i=i+1;
            if (mintraces-1)<i && i<(maxtraces+1)
                figure
                spotsTrace=spots(spots(:,10)==q,5);
                scatter(1:length(spotsTrace),spotsTrace,'s','filled')
                hold on
                line(1:length(XX)-1,XX(2:end),'lineWidth',2)
                xlim([0 max(length(XX),length(spotsTrace))])
            elseif i>maxtraces
                break
            end
        end
    catch
        disp('error')
    end
end

PwD=pdist(AllXX(1:end)');
[counts, x]=hist(PwD,1:max(PwD));
figure
histogram(PwD);
[power_spectrum_x power_spectrum_y spectrum_peaks_x spectrum_peaks_y] = FourierAndFindPeaks(x,counts,1,10000);
end