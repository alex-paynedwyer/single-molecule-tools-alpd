%select analysis sample folder, loops through cell folders

DataDir=uigetdir;
cd(DataDir)
CellDir=dir('*e*'); CellDir = CellDir(~contains({CellDir.name},'DS_Store'));
CellNum=size(CellDir);
AllSpots1=[];
AllSpots2=[];
channel=1;
SNR_thresh=0.5;
tracksonly=0;
applysegmentation=0;

for l=1:CellNum(:)
    %   CellDir(l).name
    cd(CellDir(l).name)
    %badFiles=dir('._t_*');
    %delete(badFiles.name);
    if applysegmentation==1
        matFiles=dir('*.mat');
    else
        matFiles=dir('*TRACKS.mat');
    end
    if isempty(matFiles)
        cd ..
        continue
    end
    % loop through matfiles and find one with spot data
%     for f=1:length(matFiles)
%         m = matfile(matFiles(f).name);
%         if ~isempty(whos(m,'SpotsCh1'))
%             continue
%         end
%     end
    
    for f=1:length(matFiles)
        m = matfile(matFiles(f).name);
        if ~isempty(whos(m,'SpotsCh1'))
            m1=m;
            tracksFile=matFiles(f).name;
        end
        if ~isempty(whos(m,'CellObject'))
            m2=m;
            CellObject=m2.CellObject;
        end
    end
    
    cd ..
    for ch=channel
        switch ch
            case 1
                SpotsCh1=m1.SpotsCh1;
                %SpotsCh1=SpotsCh1(SpotsCh1(:,1)>15,:); %crop LHS bright laser patches
                if isempty(SpotsCh1)
                else
                    if exist('CellObject')
                        %size(SpotsCh1)
                        SpotsCh1=SpotCellCropper(SpotsCh1,CellObject);
                        %size(SpotsCh1)
                    end
                    SpotsCh1=SpotsCh1(SpotsCh1(:,11)>SNR_thresh,:);
                    if tracksonly==1
                        SpotsCh1=SpotsCh1(SpotsCh1(:,10)>0,:); 
                    end
                    AllSpots1=cat(1,AllSpots1,SpotsCh1);
                end
            case 2
                if ~isempty(whos(m1,'SpotsCh2'))
                    SpotsCh2=m1.SpotsCh2;
                    if isempty(SpotsCh2)
                    else
                        if exist('CellObject')
                            SpotsCh2(:,1)=SpotsCh2(:,1);
                            SpotsCh2(:,1)=SpotsCh2(:,1)-size(CellObject,2);                            
                            %size(SpotsCh2)
                            SpotsCh2=SpotCellCropper(SpotsCh2,CellObject);
                            %size(SpotsCh2)
                        end                          
                        SpotsCh2=SpotsCh2(SpotsCh2(:,11)>SNR_thresh,:);    
                        if tracksonly==1
                        SpotsCh2=SpotsCh2(SpotsCh2(:,10)>0,:); 
                        end
                        AllSpots2=cat(1,AllSpots2,SpotsCh2);
                    end
                else
                end
        end
        
    end
    
    
end
for ch=channel
    switch ch
        case 1
            Spots=AllSpots1;
            %startframe = 2340;
            %endframe = 2500;
            %Spots=AllSpots1(AllSpots1(:,9)>=startframe & AllSpots1(:,9)<=endframe,:); % select range of frames to display
            
            if isempty(Spots)
                disp('nothing in channel 1')
            else
                disp(strcat("Number of channel 1 spots: ",num2str(length(AllSpots1))))
                disp(strcat("Number of channel 1 tracks: ",num2str(length(unique(AllSpots1(:,10))))))
            end
        case 2
            Spots=AllSpots2;
            if isempty(Spots)
                disp('nothing in channel 2')
            else
                disp(strcat("Number of channel 2 spots: ",num2str(length(AllSpots2))))
                disp(strcat("Number of channel 2 tracks: ",num2str(length(unique(AllSpots2(:,10))))))
            end
    end
        
    figure;
    
    subplot(2,4,1)
    [counts, x]=KDFplot(Spots(Spots(:,10)>0,5));
    xlabel('Intensity')
    ylabel('Probability Density')
    [pks,locs] = findpeaks(counts,x);
    hold on
    text(locs(pks==max(pks)), pks(pks==max(pks))*1.05,num2str(locs(pks==max(pks)),2))
    [counts2, x2]=KDFplot(Spots(Spots(:,10)>0 & Spots(:,9)>max(Spots(:,9))/3,5));
    [pks2,locs2] = findpeaks(counts2,x2);
    text(locs2(pks2==max(pks2)), pks2(pks2==max(pks2))*1.05,num2str(locs2(pks2==max(pks2)),2))
    legend('All Spots','All Spots after 1/3 bleach')
    title('Complete range')
    
    subplot(2,4,2)
    KDFplot(Spots(Spots(:,10)>0 & Spots(:,5)<(locs(pks==max(pks))*5),5));
    xlabel('Intensity')
    ylabel('Probability Density')
    hold on
    KDFplotPeaks(Spots(Spots(:,10)>0 & Spots(:,9)>max(Spots(:,9))/3 & Spots(:,5)<(locs(pks==max(pks))*5),5));
    legend('All Spots','found peaks','All Spots after 1/3 bleach')
    title('Likely Isingle range')
    
    subplot(2,4,3)
    try
        PwD=pdist(Spots(Spots(:,10)>0 & Spots(:,9)<max(Spots(:,9))/3,5));
        [counts, x]=hist(PwD,1:max(PwD));
        bar(x(x>1),counts(x>1))
        title('Pairwise distance distribution, first third of time points')
    catch
        try
            PwD=pdist(Spots(Spots(:,10)>0 & Spots(:,9)<max(Spots(:,9))/10,5));
            [counts, x]=hist(PwD,1:max(PwD));
            bar(x(x>1),counts(x>1))
            title('Pairwise distance distribution, first tenth of time points')
        catch
            PwD=pdist(Spots(Spots(:,10)>0 & Spots(:,9)<max(Spots(:,9))/100,5));
            [counts, x]=hist(PwD,1:max(PwD));
            bar(x(x>1),counts(x>1))
            title('Pairwise distance distribution, first hundredth of time points')
        end
    end
    xlabel('Intensity Step')
    ylabel('Probability Density')
    xlim([0,locs(pks==max(pks))*5])
    
    subplot(2,4,4)
    [power_spectrum_x, power_spectrum_y, spectrum_peaks_x, spectrum_peaks_y] = FourierAndFindPeaks(x,counts,0,locs(pks==max(pks))*5);
    plot(power_spectrum_x, power_spectrum_y)
    xlabel('Intensity Step')
    ylabel('Power')
    xlim([0,locs(pks==max(pks))*5])
    
    subplot(2,4,5)
    scatter(Spots(Spots(:,10)>0,9),Spots(Spots(:,10)>0,5))
    ylabel('Intensity')
    xlabel('Time (frames)')
    title('Complete range')
    
    subplot(2,4,6)
    scatter(Spots(Spots(:,10)>0,9),Spots(Spots(:,10)>0,5))
    ylabel('Intensity')
    xlabel('Time (frames)')
    ylim([0,(locs(pks==max(pks))*5)])
    hold on
    plot(0:max(Spots(Spots(:,10)>0,9)),ones(max(Spots(Spots(:,10)>0,9))+1,1)*(locs(pks==max(pks))),'--')
    title('Likely Isingle range')
       
    subplot(2,4,7)
    spotSNRhist=histogram(Spots(:,11),'BinWidth',0.01);
    xlabel('SNR');
    ylabel('Number of spots');
    xlim(gca,[0.3 1]);
        
    subplot(2,4,8)
    [xData, yData] = prepareCurveData(spotSNRhist.BinEdges(1:end-1), spotSNRhist.BinCounts);  % Set up fittype and options.
    ft = fittype( 'a*x^(-3)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Robust = 'LAR';
    opts.StartPoint = 100;
    [fitresult, gof] = fit( xData, yData, ft, opts );% Fit model to data.
    h = plot( fitresult, xData, yData );
    set(gca,'XMinorTick','on','XScale','log','YMinorTick','on','YScale','log');
    legend(gca,'off');
    xlabel('SNR');
    ylabel('Number of spots');
    xlim(gca,[0.1 10]);
    ylim(gca,[1 10000]);
    box on;
    grid on;
        
end