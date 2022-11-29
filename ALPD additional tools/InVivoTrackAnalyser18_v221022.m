%select analysis sample folder, loops through cell folders

%%Useful visualisation functions
%scatter(SpotsCh1(:,1),SpotsCh1(:,2)); hold on; scatter(SpotsCh2(:,1)-256,SpotsCh2(:,2));
%ans1=ans; figure; scatter(ans1.AllSpots1linked(:,1),ans1.AllSpots1linked(:,2)); hold on; scatter(ans1.AllSpots2linked(:,1),ans1.AllSpots2linked(:,2));

function [output,h1,h2]=InVivoTrackAnalyser16(params,h1,h2)
%clearvars -except h1 h2

%% Input parameters
if ~exist('params')  %careful - if params already exists in memory, that will be used instead of the below
    
    %SAVE/OVERWRITE OUTPUT
    params.saveOutput=1; %save to output.mat
    
    %SPOT ACCEPTANCE CRITERIA and SINGLE MOLECULE SCALING
    % note on M1, typically Ch1 = RED fluorescence (561nm), Ch2 = GREEN/YELLOW (488/514nm)
    
    params.IsingleCh1=130;%82; %133;%157; +- ~10% mCherry in vivo M1 141021 (M1 Prime95b @ 40 mW 488nm source, 5ms exposure, 7ms framerate, All Rows, 12-bit Sensitivity mode)
    params.IsingleCh2=115;%90; %114;%175; +- ~10% mGFP in vivo M1 141021 (M1 Prime95b @ 20 mW 561nm source, 5ms exposure, 7ms framerate, All Rows, 12-bit Sensitivity mode)
    params.SNRminCh1=0.2; %M1 red %default=0.4. Spots with lower SNR will be discarded. Set to the noise floor.
    params.SNRminCh2=0.2; %M1 green/yellow
    params.overtrack=0;  %extends track analysis to timepoints past SNRmin; used only for 'overtracking' to extract sm steps.
        
    %Ch1: 174-237 Ch2: 219-242 (DnaT- 75mm)
    %Ch1: 249!-274 Ch2: 219-237!-253 (DnaT+ 75mm)
    
    %select to do various extra analysis
    %=0 plots nothing!
    %=1 normal Stoichiometry and diffusion etc..
    %=2 colocalisation analysis different plots
    %=3 colocalisation analysis same plots
    %=4 scatter plot G vs R
    %=5 periodicity analysis
    %=6 CDF analysis
    %=7 hotplot
    %=8 coloc size analysis
    %=9 ??
    %=10 colocalisation liklihood
    
    params.plotSelect=3;
    %params.plotSelect=2;
    %params.plotSelect=5;
    
    %kernel density widths for plots for stoich and diffusion distributions
    params.stoichKDW = 0.6;%molecules   %0.7 standard estimated from fraction of noise on the intensity step size for a single molecule
    params.diffKDW = 0.08;%um2/s   %depends on framerate: estimate = (localisation precision)^2/(time between like frames)
    
    % Channel row width (e.g, Prime 95B is 1200x1200 so each channel is 600 wide)
    params.framewidth = 600;%600;
    
    %SPATIAL REGISTRATION of CHANNELS
    %Adjusted offsets for registering fluorescence channels. 
    % Typically 0 is fine since brightfield registration does most of the work
    xoffset = 0;
    yoffset = 0;
    %params.alignmentTransform=[xoffset-params.framewidth,-yoffset]; %for general offsets from ImageMakerAligned;
    
    %if 0, then each segment in each field analysed in turn (default)
    %if 1, then all segments in each field added together.
    %if 2, then nucleus (NuclObject) only.  (untested)
    %if 3, then all segments excluding nucleus. (untested)
    params.compartmentSwitch=0;
    
    %APPLYING SEGMENTATION
    % runs without segmentation
    params.noSeg=0;
    % warps segments generated on brightfield channel to fit same fluorescence channel
    params.warpSeg=1;
    % shifts segmentation from fluorescence channel 2 to fluorescence channel 1
    params.circshiftSeg=1; 
    % dilates segmentation with disk structural element of this size, default=2
    params.padSeg=2;
    % ad hoc window segmentation
    params.exclude_left = 0; % Parameters to exclude spots from the edges;
    params.exclude_right = 0;
    
    %COLOCALISING TRACKS and SCALING
    params.colocSwitch=1;
    %=1, links all ch1 spots to ch2
    %=2, links all ch2 spots to ch1
    
    params.lowertracklimit = 3;  % shortest/longest track lengths to include for stoichiometry analysis
    params.uppertracklimit = 15;%15; %could be 25 for GFP (if 75ms is decay time then 1/3 of molecules bleached in 15frames x2ms exposure
    
    params.frameLimit=10;%10; %number of frames to use (typically first 75ms after laser on, to avoid photobleached stoichiometry)
    %300 %= number of 10ms frames causing bleaching of GFP to 80% with 1mW laser on expanded M4 beam. 
    params.pixelSize=0.053;%0.053;%53; %pixel size in microns (11um pixel / 100x objective / 2x telescope)
    params.frameTime=0.0076;%0.0046; %frame time in seconds  (e.g. 5 ms exposure + readout of 2.5ms)
    params.PSFsize=200; %size of PSF in nm
    % any co-ordinate transformation between the two channels
    params.transform=[0,0];
    params.colocError=2.5; %max colocalisation error in pixels ~3x localisation precision
    
    % maximum distance in pixels to link spots across channels
    params.d=7;
    % minumum overlap integral to link spots across channels
    params.overlap=0.75; %0.75 is highest sensible value. 0.2 is lowest possible = Rayleigh criterion; 0.9 is far too high.
        
    % set if=0 link all spots regardless of frame
    %     if=1 link only spots in the same frame
    %     if=2 link spots in alternate frames (for ALEX)
    params.frameLinkMethod=2;
    
    % set so that spots are only assigned 1 partner
    params.Unique=0;
    
    %PER CELL STATISTICS
    %set to display cell specific colocalisation graphs
    params.showOutput=0;
    
    %set to 1 to ignore cells empty in both channels, currently just cuts
    %out of the means number of tracks/cell
    params.IgnoreEmpty=1;
    
    try
        params.DataDir=uigetdir;
        cd(params.DataDir)
        DataDir=params.DataDir
    catch
        error("No DataDir directory selected.")
    end
end

CellDir=dir('*ideo*'); %Sonam's convention 'video'
%CellDir=dir('*el*'); %Standard for 'cell' and 'field'
CellNum=size(CellDir);
AllSpots1=[];
AllSpots2=[];
AllSpots1linked=[];
AllSpots2linked=[];
trackArrayCh1=[];
trackArrayCh2=[];
spotsInTracksCh1=[];
spotsInTracksCh2=[];
segArray=[];
fieldNo=0;
cellNo1=0;
cellNo2=0;
emptyCh1=1;
emptyCh2=1;

if params.overtrack==1
    matFiles=dir('*BASELINE.mat');
elseif params.overtrack==0
    matFiles=dir('*TRACKS.mat');
end

if isempty(matFiles)
    [sub,fls] =subdir;
else
    fls={1};
    sub={pwd};
end

for l=1:length(fls)
    if ~isempty(fls{l})
        cd(sub{l})
        if params.overtrack==1
            matFiles=[dir('*BASELINE.mat'),dir('*segmentation.mat')];
        elseif params.overtrack==0
            matFiles=[dir('*TRACKS.mat'),dir('*segmentation.mat')];
        end
        % loop through matfiles and find one with spot data
        for f=1:length(matFiles)
            m = matfile(matFiles(f).name);
            if ~isempty(whos(m,'SpotsCh1'))
                m1=m;
                tracksFile=matFiles(f).name;
            end
            if ~isempty(whos(m,'CellObject')) % || exist(m,'CellObject')
                m2=m;
                fieldNo=fieldNo+1;
                %  segObject=m2.CellObject;
                if params.compartmentSwitch==0
                    segObject=m2.CellObject;
                elseif params.compartmentSwitch==1
                    CellObject=sum(m2.CellObject,3);
                    CellObject(CellObject>1)=1;
                    segObject=CellObject;
                elseif params.compartmentSwitch==2
                    segObject=NucObject;
                elseif params.compartmentSwitch==3
                    NucObject=sum(m2.NucObject,3);
                    NucObject(NucObject>1)=1;
                    segObject=CellObject-NucObject;
                end
                segObject(segObject<0)=0;
                if params.padSeg>0
                    segObject=imdilate(segObject, strel('disk',params.padSeg));
                end
                
            end
        end
        
        if ~exist('segObject') || isempty(segObject)
            % disp('fake seg')
            if params.noSeg==1
                %spots=m1.SpotsCh2;
                
                %EDIT 1
                %segObject=ones(260,375);
                segObject=ones(size(m1.frame_average));  %could cause issues if segmentation area is smaller than actual image (e.g. empty channel cropped out)
                
                %segObject=ones(ceil(max(spots(:,2))),ceil(max(spots(:,1))));
            else
                continue
            end
        end
        
        %  cd ..
        for ch=1:2
            switch ch
                case 1
                    if ~isempty(whos(m1,'SpotsCh1')) && length(m1.SpotsCh1)>=12
                    SpotsCh1=m1.SpotsCh1;
                    %EDIT 3A
                    segObjectCh1=segObject;
                    tform=m1.tform;
                        if params.warpSeg==1
                        segObjectCh1 = imwarp(segObjectCh1,tform,'OutputView',imref2d(size(segObjectCh1)));
                        end
                        if params.circshiftSeg==1
                        segObjectCh1 = circshift(segObjectCh1,[0,size(segObjectCh1,2)/2]);
                        end
                    % assign image number
                    if isempty(SpotsCh1) || size(SpotsCh1,2)<12
                        cellNo1=cellNo1+size(segObjectCh1,3);
                        disp("No spots in Channel 1.")
                        continue
                    else
                        SpotsCh1=SpotsCh1(SpotsCh1(:,1)>params.exclude_left & SpotsCh1(:,1)<(params.framewidth-params.exclude_right),:); %crop edges
                        %EDIT 2A
                        SpotsCh1=SpotsCh1(SpotsCh1(:,11)>params.SNRminCh1,:); %screen low SNR spots
                        
                        
                        if SpotsCh1(1,5)==0
                            cellNo1=cellNo1+size(segObject,3);
                            continue
                        else
                            spots=SpotsCh1;
                            SpotsCh1temp =spots;
                            SpotsCh1temp(:,13)=l;
                            AllSpots1=cat(1,AllSpots1,SpotsCh1temp);
                            [trackArrayCh1,spotsInTracksCh1,segArray,cellNo1]=trackAnalyser5_ALH_general(spots,segObjectCh1,trackArrayCh1,params.frameLimit,params.frameTime,params.pixelSize,tracksFile,params.PSFsize,segArray,cellNo1,params.lowertracklimit,params.uppertracklimit);
                        end
                        emptyCh1=0;
                    end
                    end
                case 2
                    if ~isempty(whos(m1,'SpotsCh2')) && length(m1.SpotsCh2)>=12
                        SpotsCh2=m1.SpotsCh2;
                        %EDIT 3B
                        tform=m1.tform;
                        if params.warpSeg==1
                        segObjectCh2 = imwarp(segObject,tform,'OutputView',imref2d(size(segObject)));
                        end
                        % assign image number
                        SpotsCh2=SpotsCh2(SpotsCh2(:,1)>params.framewidth+params.exclude_left & SpotsCh2(:,1)<(2*params.framewidth-params.exclude_right-1),:);
                        %EDIT 2B
                        SpotsCh2=SpotsCh2(SpotsCh2(:,11)>params.SNRminCh2,:); %screen low SNR spots
                        if isempty(SpotsCh2)  || size(SpotsCh2,2)<12
                            
                            cellNo2=cellNo1;
                            continue
                        else
                            if SpotsCh2(1,5)==0
                                cellNo2=cellNo1;
                                continue
                            else
                                %AllSpots2=cat(1,AllSpots2,SpotsCh2temp);
                                spots=SpotsCh2;
                                params.alignmentTransform=[-xoffset-params.framewidth,-yoffset];
                                spots(:,1)=spots(:,1)+params.alignmentTransform(1);
                                spots(:,2)=spots(:,2)+params.alignmentTransform(2);

                                SpotsCh2temp =spots;
                                SpotsCh2temp(:,13)=l;
                                AllSpots2=cat(1,AllSpots2,SpotsCh2temp);
                                segObjectCh2=circshift(segObject,[-params.alignmentTransform(2),-params.alignmentTransform(1)]);
                                %params.segObjectCh2=segObjectCh2;
                                [trackArrayCh2,spotsInTracksCh2,~,cellNo2]=trackAnalyser5_ALH_general(spots,segObjectCh2,trackArrayCh2,params.frameLimit,params.frameTime,params.pixelSize,tracksFile,params.PSFsize,segArray,cellNo2,params.lowertracklimit,params.uppertracklimit);
                            end
                        end
                        emptyCh2=0;
                    else
                        cellNo2=cellNo1;
                        disp("No spots in Channel 2.")
                        continue
                    end
            end
        end
        if ~isempty(spotsInTracksCh1) && ~isempty(spotsInTracksCh2)
            %   if ~isempty(SpotsCh1(SpotsCh1(:,10)>0,:)) && ~isempty(SpotsCh2(SpotsCh2(:,10)>0,:)) && SpotsCh1(1,5)>0  && SpotsCh2(1,5)>0
            %    try
            if params.colocSwitch==1
                [SpotsCh1linked, SpotsCh2linked]=Colocaliser2(spotsInTracksCh1,spotsInTracksCh2,params);
            else
                [SpotsCh2linked,SpotsCh1linked]=Colocaliser2(spotsInTracksCh2,spotsInTracksCh1,params);
            end
            % this is section is a loop for now, possible with arrays but
            % difficult, safer to stick to loop for now
            [ch1Traj,ch1Ind]=unique(SpotsCh1linked(SpotsCh1linked(:,14)>0,10)); % linked trajectories in ch1
            ch2Traj=SpotsCh1linked(SpotsCh1linked(:,14)>0,15);
            for tr=1:length(ch1Traj) %loop over linked trajectory numbers
                % link tractory in ch1 with the longest linked traj in ch2
                ch2Traj(tr)=mode(SpotsCh1linked(SpotsCh1linked(:,14)>0 & SpotsCh1linked(:,10)==ch1Traj(tr),15));
                trackArrayCh1(trackArrayCh1(:,7)==trackArrayCh1(end,7) & trackArrayCh1(:,6)==ch1Traj(tr),8)=ch2Traj(tr);
                % length of time in frames linked to a trajectory
                trackArrayCh1(trackArrayCh1(:,7)==trackArrayCh1(end,7) & trackArrayCh1(:,6)==ch1Traj(tr),9)=sum(SpotsCh1linked(:,10)==ch1Traj(tr) & SpotsCh1linked(:,15)==ch2Traj(tr));
                try
                    trackArrayCh1(trackArrayCh1(:,7)==trackArrayCh1(end,7) & trackArrayCh1(:,6)==ch1Traj(tr),10)=find(trackArrayCh2(:,7)==trackArrayCh2(end,7) & trackArrayCh2(:,6)==ch2Traj(tr));
                catch
                    trackLength = find(trackArrayCh2(:,7)==trackArrayCh2(end,7) & trackArrayCh2(:,6)==ch2Traj(tr));
                    trackArrayCh1(trackArrayCh1(:,7)==trackArrayCh1(end,7) & trackArrayCh1(:,6)==ch1Traj(tr),10)= trackLength(1);
                end
                % distance from other traj
                trackArrayCh1(trackArrayCh1(:,7)==trackArrayCh1(end,7) & trackArrayCh1(:,6)==ch1Traj(tr),11)=mean(SpotsCh1linked(SpotsCh1linked(:,10)==ch1Traj(tr) & SpotsCh1linked(:,15)==ch2Traj(tr),16));
                % link corresponding tractory in ch2 with ch1
                trackArrayCh2(trackArrayCh2(:,7)==trackArrayCh2(end,7) & trackArrayCh2(:,6)==ch2Traj(tr),8)=ch1Traj(tr);
                % assign the same link time to ch2
                trackArrayCh2(trackArrayCh2(:,7)==trackArrayCh2(end,7) & trackArrayCh2(:,6)==ch2Traj(tr),9)=sum(SpotsCh1linked(:,10)==ch1Traj(tr) & SpotsCh1linked(:,15)==ch2Traj(tr));
                try
                    trackArrayCh2(trackArrayCh2(:,7)==trackArrayCh2(end,7) & trackArrayCh2(:,6)==ch2Traj(tr),10)=find(trackArrayCh1(:,7)==trackArrayCh1(end,7) & trackArrayCh1(:,6)==ch1Traj(tr));
                catch
                    trackLink = find(trackArrayCh1(:,7)==trackArrayCh1(end,7) & trackArrayCh1(:,6)==ch1Traj(tr));
                    trackArrayCh2(trackArrayCh2(:,7)==trackArrayCh2(end,7) & trackArrayCh2(:,6)==ch2Traj(tr),10) = trackLink(1);
                end
                trackArrayCh2(trackArrayCh2(:,7)==trackArrayCh2(end,7) & trackArrayCh2(:,6)==ch2Traj(tr),11)=mean(SpotsCh1linked(SpotsCh1linked(:,10)==ch1Traj(tr) & SpotsCh1linked(:,15)==ch2Traj(tr),16));
            end
            %             SpotsCh1linked(:,17)=max(trackArrayCh1(1));
            %             SpotsCh2linked(:,17)=max(trackArrayCh1(1));
            AllSpots1linked=cat(1,AllSpots1linked,SpotsCh1linked);
            AllSpots2linked=cat(1,AllSpots2linked,SpotsCh2linked);
            %             catch err
            %             end
            %  end
        else
            %             trackArrayCh1(trackArrayCh1(:,7)==trackArrayCh1(end,7),8:9)=0;
            %             trackArrayCh2(trackArrayCh2(:,7)==trackArrayCh2(end,7),8:9)=0;
            
        end
    end
    
    
end

%% Generate some outputs
output.AllSpots1=AllSpots1;
output.AllSpots2=AllSpots2;
output.trackArrayCh1=trackArrayCh1;
output.trackArrayCh2=trackArrayCh2;
output.AllSpots1linked=AllSpots1linked;
output.AllSpots2linked=AllSpots2linked;
output.segArray=segArray;
output.NFields=fieldNo;

%% Start the plotting process
cmap=colormap('jet');
params.cmap=cmap;
maxColor=5;
params.maxColor=maxColor;
close(gcf)
for ch=1:2
    switch ch
        case 1
            if isempty(trackArrayCh1)
                h1=[];
                continue
            end
            plotName='Track Plot Ch1';
            if exist('h1')
                figure(h1)
                hchild1=get(h1,'children');
                axisProps=findobj(hchild1(end),'Type','line');
                colorInd=size(axisProps,1)+1;
            else
                colorInd=1;
                if params.plotSelect==0
                    h1=[];
                else
                    h1=figure('Name',plotName,'NumberTitle','off');
                end
            end
            trackArray=trackArrayCh1;
            Isingle=params.IsingleCh1;
        case 2
            if isempty(trackArrayCh2)
                h2=[];
                continue
            end
            plotName='Track Plot Ch2';
            if exist('h2')
                figure(h2)
                hchild2=get(h2,'children');
                axisProps=findobj(hchild2(end),'Type','line');
                colorInd=size(axisProps,1)+1;
            else
                colorInd=1;
                if params.plotSelect==0
                    h2=[];
                else
                    h2=figure('Name',plotName,'NumberTitle','off');
                end
            end
            trackArray=trackArrayCh2;
            Isingle=params.IsingleCh2;
    end
    %% means
    disp(strcat('N fields = ',num2str(fieldNo)))
    
    if  params.IgnoreEmpty==0
        for c=1:max(trackArray(:,1))
            unlinkTrack(c)=sum(trackArray(trackArray(:,1)==c,8)==0);
            linkTrack(c)=sum(trackArray(trackArray(:,1)==c,8)>0);
            
        end
        disp(strcat('N all cells = ',num2str(max(trackArray(:,1)))))
        output.NCells(ch)=max(trackArray(:,1));
    else
        try
            uniqueCellNo=unique([trackArrayCh1(:,1)',trackArrayCh2(:,1)']);
        catch
            try
            uniqueCellNo=unique([trackArrayCh1(:,1)']);
            catch
            uniqueCellNo=unique([trackArrayCh2(:,1)']);
            end
        end
        for c=1:length(uniqueCellNo)
            unlinkTrack(c)=sum(trackArray(trackArray(:,1)==uniqueCellNo(c),8)==0);
            if unlinkTrack(c)>0
                unlinkTotal(c)=sum(trackArray(trackArray(trackArray(:,1)==uniqueCellNo(c),8)==0,2))/Isingle;
            else
                unlinkTotal(c)=0;
            end
            linkTrack(c)=sum(trackArray(trackArray(:,1)==uniqueCellNo(c),8)>0);
            if linkTrack(c)>0
                linkTotal(c)=sum(trackArray(trackArray(trackArray(:,1)==uniqueCellNo(c),8)>0,2))/Isingle;
            else
                linkTotal(c)=0;
            end
        end
        NumberSegments=length(uniqueCellNo);
        disp(strcat('N non-empty cell segments = ',num2str(NumberSegments)))
        output.NCells(ch)=length(uniqueCellNo);
    end
    
    NlinkTrack = sum(trackArray(:,8)>0);
    NunlinkTrack =sum(trackArray(:,8)==0);
    
    output.linkTrackMean(1,ch)=mean(linkTrack);
    output.linkTrackMean(2,ch)=std(linkTrack)/NumberSegments^0.5;
    output.unlinkTrackMean(1,ch)=mean(unlinkTrack);
    output.unlinkTrackMean(2,ch)=std(unlinkTrack)/NumberSegments^0.5;
    
    output.alllinkTrack(:,ch)=linkTrack;
    output.allunlinkTrack(:,ch)=unlinkTrack;
    
    output.NlinkTrack(ch)= NlinkTrack;
    output.NunlinkTrack(ch)= NunlinkTrack;
    
    output.linkStoichsMean(1,ch)=mean(trackArray(trackArray(:,8)>0,2)/Isingle);
    output.linkStoichsMean(2,ch)=std(trackArray(trackArray(:,8)>0,2)/Isingle)/(NlinkTrack^0.5);
    output.unlinkStoichsMean(1,ch)=mean(trackArray(trackArray(:,8)==0,2)/Isingle);
    output.unlinkStoichsMean(2,ch)=std(trackArray(trackArray(:,8)==0,2)/Isingle)/(NunlinkTrack^0.5);
    
    output.linkDiffsMean(1,ch)=mean(trackArray(trackArray(:,8)>0,3));
    output.linkDiffsMean(2,ch)=std(trackArray(trackArray(:,8)>0,3))/(NlinkTrack^0.5);
    output.unlinkDiffsMean(1,ch)=mean(trackArray(trackArray(:,8)==0,3));
    output.unlinkDiffsMean(2,ch)=std(trackArray(trackArray(:,8)==0,3))/(NunlinkTrack^0.5);
    output.separationMean(1,ch)=mean(trackArray(trackArray(:,11)>0,11)*params.pixelSize*1000);
    output.separationMean(2,ch)=std(trackArray(trackArray(:,11)>0,11)*params.pixelSize*1000)/(NlinkTrack)^0.5;
    
    disp(strcat('channel = ',num2str(ch)))
    
    disp(strcat('N total trajectories/cell = ',num2str((mean(unlinkTrack+linkTrack)),3),'+/-',num2str((std(unlinkTrack+linkTrack)/(NumberSegments)^0.5),3)))
    disp(strcat('N unlinked trajectories/cell = ',num2str(mean(unlinkTrack),3),'+/-',num2str(std(unlinkTrack)/(NumberSegments)^0.5,3)))
    disp(strcat('N linked trajectories/cell = ',num2str(mean(linkTrack),3),'+/-',num2str(std(linkTrack)/(NumberSegments)^0.5,3)))
    
    disp(strcat('N molecules in total trajectories/cell = ',num2str(sum(trackArray(:,2)/Isingle)/NumberSegments,3),'+/-',num2str(std(unlinkTotal+linkTotal)/(NumberSegments)^0.5,3)))
    disp(strcat('N molecules in unlinked trajectories/cell = ',num2str(sum(trackArray(trackArray(:,8)==0,2)/Isingle)/NumberSegments,3),'+/-',num2str(std(unlinkTotal)/NumberSegments^0.5,3)))
    disp(strcat('N molecules in linked trajectories/cell = ',num2str(sum(trackArray(trackArray(:,8)>0,2)/Isingle)/NumberSegments,3),'+/-',num2str(std(linkTotal)/NumberSegments^0.5,3)))
    
    disp(strcat('Mean stoichiometry = ',num2str(mean(trackArray(:,2)/Isingle),3),'+/-',num2str(std(trackArray(:,2)/Isingle)/(NunlinkTrack+NlinkTrack)^0.5,3)))
    disp(strcat('Mean unlinked stoichiometry = ',num2str(mean(trackArray(trackArray(:,8)==0,2)/Isingle),3),'+/-',num2str(std(trackArray(trackArray(:,8)==0,2)/Isingle)/(NunlinkTrack)^0.5,3)))
    disp(strcat('Mean linked stoichiometry = ',num2str(mean(trackArray(trackArray(:,8)>0,2)/Isingle),3),'+/-',num2str(std(trackArray(trackArray(:,8)>0,2)/Isingle)/(NlinkTrack)^0.5,3)))
    
    disp(strcat('Mean D = ',num2str(mean(trackArray(:,3)),2),'+/-',num2str(std(trackArray(:,3))/(NunlinkTrack+NlinkTrack)^0.5,2)))
    disp(strcat('Mean unlinked D (um2/s) = ',num2str(mean(trackArray(trackArray(:,8)==0,3)),2),'+/-',num2str(std(trackArray(trackArray(:,8)==0,3))/(NunlinkTrack)^0.5,2)))
    disp(strcat('Mean linked D (um2/s) = ',num2str(mean(trackArray(trackArray(:,8)>0,3)),2),'+/-',num2str(std(trackArray(trackArray(:,8)>0,3))/(NlinkTrack)^0.5,2)))
    %disp(strcat('Mean linked separation = ',num2str(mean(trackArray(trackArray(:,11)>0,11)*params.pixelSize*1000),3),'+/-',num2str(std(trackArray(trackArray(:,11)>0,11)*params.pixelSize*1000),2)))
    
    %%%
    %% ALH extra calculations to store complete distributions outside trackArrays
    % allows multiple datasets to be collated after trajectory analysis, even with different params (Isingles, framerates)
    % simply by concatenating the relevant lists from the output files
    
    output.CellNoofTrajList(ch) = {trackArray(:,1)};
    
    stoichs = trackArray(:,2)/Isingle;
    unlinkedstoichs = trackArray(trackArray(:,8)==0,2)/Isingle;
    linkedstoichs = trackArray(trackArray(:,8)>0,2)/Isingle;
    
    if emptyCh1+emptyCh2==0 && sum(linkedstoichs) > 0
    pairedstoichs(:,1) = trackArrayCh1(trackArrayCh1(:,8)>0,2)/params.IsingleCh1;
    pairedstoichs(:,2) = trackArrayCh2(trackArrayCh1(trackArrayCh1(:,8)>0,10),2)/params.IsingleCh2;
    else
    pairedstoichs = [];
    end
    
    output.UnlinkedStoichsList(ch)={unlinkedstoichs};
    output.LinkedStoichsList(ch)={linkedstoichs};
    output.PairedStoichsList = pairedstoichs;
    
    diffs = trackArray(:,3);
    unlinkeddiffs=trackArray(trackArray(:,8)==0,3);
    linkeddiffs=trackArray(trackArray(:,8)>0,3);
    if emptyCh1+emptyCh2==0 && sum(linkedstoichs) > 0
    paireddiffs(:,1) = trackArrayCh1(trackArrayCh1(:,8)>0,3);
    paireddiffs(:,2) = trackArrayCh2(trackArrayCh1(trackArrayCh1(:,8)>0,10),3);
    else
    paireddiffs = [];
    end
    
    output.UnlinkedDiffsList(:,ch)={unlinkeddiffs};
    output.LinkedDiffsList(:,ch)={linkeddiffs};
    output.PairedDiffsList = paireddiffs;
    
    unlinkedtrajlist = trackArray(trackArray(:,8)==0,1);
    linkedtrajlist = trackArray(trackArray(:,8)>0,1);
    
    output.UnlinkedTrajList(:,ch)={unlinkedtrajlist};
    output.LinkedTrajList(:,ch)={linkedtrajlist};
    
    colourmap = cmap(round(64*(colorInd)/maxColor),:);
    
    %% Hardcore Analysis
    switch params.plotSelect
        case 0
            
        case 1
            % plots
            
            subplot(2,4,1)
            %figure('Name','Stoichiometry Distribution','NumberTitle','off')
            [plotH,~,x]=KDFplotH(stoichs,params.stoichKDW);
            set(plotH,'Color',colourmap);
            xlabel('Stoichiometry')
            ylabel('Probability Density')
            hold on
            xlim([0,max(x)])
            pbaspect([1 1 1]);
            
            subplot(2,4,2)
            %figure('Name','Stoichiometry Histogram','NumberTitle','off')
            histogram(stoichs,0.5:1:round(max(stoichs)+0.5));%'FaceColor',colourmap);
            xlabel('Stoichiometry')
            ylabel('Frequency')
            xlim([0,max(stoichs)*1.1])
            pbaspect([1 1 1]);
            hold on
            
            subplot(2,4,3)
            %figure('Name','Diffusion Const Distribution','NumberTitle','off')
            [plotH,~,~]=KDFplotH(diffs,params.diffKDW);
            set(plotH,'Color',colourmap)
            pbaspect([1 1 1]);
            xlabel('Diffusion Coefficient ({\mu}m^{2}/s)')
            ylabel('Probability Density')
            xlim([0,max(diffs)*1.1])
            hold on
            
            subplot(2,4,4)
            %figure('Name','S vs D','NumberTitle','off')
            %scatter(AllDiff, AllmeanI/Isingle);
            scatter(diffs,stoichs,50,colourmap,'filled');
            xlabel('Diffusion Coefficient ({\mu}m^{2}/s)')
            ylabel('Stoichiometry')
            box on
            xlim([0,max(diffs)*1.1])
            ylim([0,max(stoichs)*1.1])
            pbaspect([1 1 1]);
            hold on
            
            %no spots/cell
            subplot(2,4,5)
            [counts,segedges]=histcounts(trackArray(:,1),'BinWidth',1);
            histogram(counts,'FaceColor',colourmap);
            xlabel('Number of trajectories/cell')
            ylabel('Frequency')
            xl = xlim;
            xlim([-0.5,xl(2)])
            pbaspect([1 1 1]);
            hold on
            
            %total number of molecules in spots/segment
            subplot(2,4,6)
            for c=1:length(trackArray(:,1))
                molSpotPerCell(c)=sum(trackArray(trackArray(:,1)==c,2)/Isingle);
            end
            [plotH,~,x]=KDFplotH(molSpotPerCell,params.stoichKDW);
            set(plotH,'Color',colourmap);
            xlim([0,max(x)])
            xlabel('total number of molecules in trajectories/cell')
            ylabel('Probability Density')
            pbaspect([1 1 1]);
            hold on
            
            %length of tracks
            subplot(2,4,7)
            for c=1:length(trackArray(:,1))
                tracklengths(c)=trackArray(c,14);
                tracktimes(c)=tracklengths(c)*params.frameTime*1000;
            end
            [plotH,~,x]=KDFplotH(tracktimes);
            xlim([0,max(x)])
            %set(plotH,'Color',colourmap);
            xlabel('Track length (ms)')
            ylabel('Probability Density')
            pbaspect([1 1 1]);
            hold on
            
            %mean spot stoich/cell
            %subplot(2,4,7)
            %for c=1:max(trackArray(:,1))
            %    molSpotPerCell(c)=mean(trackArray(trackArray(:,1)==c,2)/Isingle);
            %end
            %[plotH,~,x]=KDFplotH(molSpotPerCell);
            %xlim([0,max(x)])
            %set(plotH,'Color',colourmap);
            %xlabel('Mean number of molecules in trajectory/cell')
            %ylabel('Probability Density')
            %pbaspect([1 1 1]);
            %hold on
            
        case 2
            if min(trackArray(:,8)==0)
                disp('No linked spots found.')
            else
                %% plots
                
                subplot(2,4,1)
                %figure('Name','Stoichiometry Distribution','NumberTitle','off')
                [plotH,~,x]=KDFplotH(unlinkedstoichs,params.stoichKDW);
                set(plotH,'Color','blue')
                xlabel('Stoichiometry')
                ylabel('Probability Density')
                title('Unlinked Trajectories')
                xlim([0,max(x)])
                pbaspect([1 1 1]);
                hold on
                
                subplot(2,4,2)
                %figure('Name','Stoichiometry Distribution','NumberTitle','off')
                [plotH,~,x]=KDFplotH(linkedstoichs,params.stoichKDW);
                set(plotH,'Color','red')
                xlabel('Stoichiometry')
                ylabel('Probability Density')
                title('Linked Trajectories')
                xlim([0,max(x)])
                hold on
                pbaspect([1 1 1]);
                
                subplot(2,4,3)
                %figure('Name','Diffusion Const Distribution','NumberTitle','off')
                [plotH,~,~]=KDFplotH(unlinkeddiffs,params.diffKDW);
                set(plotH,'Color','blue')
                xlabel('Diffusion Coefficient ({\mu}m^{2}/s)')
                ylabel('Probability Density')
                title('Unlinked Trajectories')
                pbaspect([1 1 1]);
                xlim([0,max(diffs)*1.2])
                hold on
                
                subplot(2,4,4)
                %figure('Name','Diffusion Const Distribution','NumberTitle','off')
                [plotH,~,x]=KDFplotH(linkeddiffs,params.diffKDW);
                set(plotH,'Color','red')
                xlabel('Diffusion Coefficient ({\mu}m^{2}/s)')
                ylabel('Probability Density')
                title('Linked Trajectories')
                pbaspect([1 1 1]);
                xlim([0,max(diffs)*1.2])
                hold on
                
                subplot(2,4,5)
                [unlinkcounts,unlinksegedges]=histcounts(unlinkedtrajlist,'BinWidth',1);
                histogram(unlinkcounts,'FaceColor','blue')
                xlabel('Number of unlinked trajectories/cell')
                ylabel('Frequency')
                pbaspect([1 1 1]);
                hold on
                
                subplot(2,4,6)
                [linkcounts,linksegedges]=histcounts(linkedtrajlist,'BinWidth',1);
                histogram(linkcounts,'FaceColor','red')
                xlabel('Number of linked trajectories/cell')
                ylabel('Frequency')
                pbaspect([1 1 1]);
                hold on
                
                subplot(2,4,7)
                [plotH,~,x]=KDFplotH(trackArray(trackArray(:,11)>0,11)*params.pixelSize*1000,10);
                set(plotH,'Color','red')
                xlabel('separation between linked spots (nm)')
                ylabel('Probability Density')
                title('Linked Trajectories')
                pbaspect([1 1 1]);
                xlim([0,max(x)])
                hold on
                
                subplot(2,4,8)
                title('Linked Trajectories')
                scatter(trackArray(trackArray(:,11)>0,2)/Isingle,trackArray(trackArray(:,11)>0,11)*params.pixelSize*1000,50,'red','filled');
                ylabel('Separation between linked spots (nm)')
                xlabel('Stoichiometry');
                pbaspect([1 1 1]);
                hold on
            end
        case 3
            % plots
            if min(trackArray(:,8)==0)
                disp('No linked spots found.')
            else
                
                subplot(2,2,1);
                box on
                %figure('Name','Stoichiometry Distribution','NumberTitle','off')
                [plotH,~,x]=KDFplotH(unlinkedstoichs,params.stoichKDW);
                set(plotH,'Color','blue')
                xlabel('Stoichiometry')
                ylabel('Probability Density')
                %title('Unlinked Trajectories')
                hold on
                [plotH,~,x]=KDFplotH(linkedstoichs,params.stoichKDW);
                set(plotH,'Color','red')
                xlabel('Stoichiometry')
                ylabel('Probability Density')
                %title('Linked Trajectories')
                hold on
                axis tight
                pbaspect([1 1 1]);
                xlim([0,max(x)])
                
                subplot(2,2,2)
                box on
                [plotH,~,x]=KDFplotH(unlinkeddiffs,params.diffKDW);
                set(plotH,'Color','blue')
                hold on
                [plotH,~,x]=KDFplotH(linkeddiffs,params.diffKDW);
                set(plotH,'Color','red');
                xlabel('Diffusion Coefficient ({\mu}m^{2}/s)')
                ylabel('Probability Density')
                %title('linked Trajectories')
                axis tight
                pbaspect([1 1 1]);
                xlim([0,max(diffs*1.07)])
                hold on
                
                subplot(2,2,3)
                box on
                [unlinkcounts,unlinksegedges]=histcounts(unlinkedtrajlist,'BinWidth',1);
                histogram(unlinkcounts,'FaceColor','blue');
                hold on
                [linkcounts,linksegedges]=histcounts(linkedtrajlist,'BinWidth',1);
                histogram(linkcounts,'FaceColor','red');
                xlabel('Number of trajectories/cell')
                ylabel('Frequency')
                axis tight
                pbaspect([1 1 1]);
                hold on
                
                subplot(2,2,4)
                scatter(pairedstoichs(:,1),pairedstoichs(:,2),'filled')
                xlabel('Ch1 linked stoichiometry')
                ylabel('Ch2 linked stoichiometry')
                gradient=pairedstoichs(:,1)\pairedstoichs(:,2)
                hold on
                plot([0;pairedstoichs(:,1);max(pairedstoichs(:,1))*1.1],[0;pairedstoichs(:,1)*gradient;max(pairedstoichs(:,1))*1.1*gradient],'--','Color','black')
                disp(strcat('Gradient = ',num2str(gradient)))
                axis tight
                box on
                xlim([0,max(pairedstoichs(:,1))*1.1])
                ylim([0,max(pairedstoichs(:,2))*1.1])
                pbaspect([1 1 1]);
                text(max(pairedstoichs(:,1))*1,max(pairedstoichs(:,2))*1,strcat('Gradient: ',num2str(round(gradient,2))),'HorizontalAlignment','right');
                
            end
            
        case 4
            scatter(pairedstoichs(:,1),pairedstoichs(:,2),'filled')
            xlabel('Ch1 linked stoichiometry')
            ylabel('Ch2 linked stoichiometry')
            gradient=pairedstoichs(:,1)\pairedstoichs(:,2);
            hold on
            plot([0;pairedstoichs(:,1);max(pairedstoichs(:,1))*1.1],[0;pairedstoichs(:,1)*gradient;max(pairedstoichs(:,1))*1.1*gradient],'--','Color','black')
            disp(strcat('Gradient = ',num2str(gradient)))
            axis tight
            box on
            xlim([0,max(pairedstoichs(:,1))*1.1])
            ylim([0,max(pairedstoichs(:,2))*1.1])
            pbaspect([1 1 1]);
            text(max(pairedstoichs(:,1))*1,max(pairedstoichs(:,2))*1,strcat('Gradient: ',num2str(round(gradient,2))),'HorizontalAlignment','right');
            
        case 5
            figure('Name','Periodicity Analysis','NumberTitle','off');
            
            subplot(2,3,1)
            [counts,x,plotH]=KDFplotPeaks(stoichs,params.stoichKDW);
            set(plotH,'Color',colourmap)
            xlabel('Stoichiometry')
            ylabel('Probability Density')
            xlim([0,max(x)])
            title('Stoichiometry Distribution')
            
            subplot(2,3,2)
            [counts,x,plotH]=KDFplotPeaks(trackArray(trackArray(:,2)<(20*Isingle),2)/Isingle,params.stoichKDW);
            set(plotH,'Color',colourmap)
            xlabel('Stoichiometry')
            ylabel('Probability Density')
            xlim([0,20])
            title('Stoichiometry Distribution up to 20')
            
            
            subplot(2,3,3)
            PwD=pdist(stoichs);
            [counts, x]=hist(PwD,0.5:0.5:max(PwD));
            bar(x,counts)
            title('Pairwise distance distribution of stoichiometries')
            xlabel('Molecule Step')
            ylabel('Probability Density')
            %xlim([0,locs(pks==max(pks))*5])
            
            subplot(2,3,5)
            [power_spectrum_x power_spectrum_y spectrum_peaks_x spectrum_peaks_y] = FourierAndFindPeaks(x,counts,0);
            plot(power_spectrum_x, power_spectrum_y)
            xlabel('Molecule Step')
            ylabel('Power')
            xlim([0,20])
            title('Fourier Spectrum Stoichiometry up to 20')
            %xlim([0,locs(pks==max(pks))*5])
            subplot(2,3,4)
            plot(power_spectrum_x, power_spectrum_y)
            xlabel('Molecule Step')
            ylabel('Power')
            xlim([0,max(x)])
            title('Fourier Spectrum Stoichiometry')
            subplot(2,3,6)
            plot(power_spectrum_x, power_spectrum_y)
            xlabel('Molecule Step')
            ylabel('Power')
            xlim([0,10])
            title('Fourier Spectrum Stoichiometry up to 10')
        case 6
            [fitresult1,fitresult2]=CDFmobility2(trackArray(:,4));
        case 7
            [heatmapdata,PlotRange]=gausshotplot(trackArrayCh1(trackArrayCh1(:,8)>0,2)/params.IsingleCh1,trackArrayCh2(trackArrayCh1(trackArrayCh1(:,8)>0,10),2)/params.IsingleCh2,1000,25,[0,30;0,30]);
        case 8
            %%
            disp(strcat('Mean unlinked Spot Size Cell Parallel (nm)=',num2str(mean(trackArray(trackArray(:,8)==0,12)),3),'+/-',num2str(std(trackArray(trackArray(:,8)==0,12)),3)))
            disp(strcat('Mean linked Spot Size Cell Parallel (nm)=',num2str(mean(trackArray(trackArray(:,8)>0,12)),3),'+/-',num2str(std(trackArray(trackArray(:,8)>0,12)),3)))
            disp(strcat('Mean unlinked Spot Size Cell Perpendicular (nm)=',num2str(mean(trackArray(trackArray(:,8)==0,13)),3),'+/-',num2str(std(trackArray(trackArray(:,8)==0,13)),3)))
            disp(strcat('Mean linked Spot Size Cell Perpendicular (nm)=',num2str(mean(trackArray(trackArray(:,8)>0,13)),3),'+/-',num2str(std(trackArray(trackArray(:,8)>0,13)),3)))
            
            subplot(2,4,1)
            %figure('Name','Stoichiometry Distribution','NumberTitle','off')
            [plotH,~,x]=KDFplotH(trackArray(trackArray(:,8)==0,12),10);
            set(plotH,'Color',colourmap)
            xlabel('Spot Size Cell Parallel (nm)')
            ylabel('Probability Density')
            title('Unlinked Trajectories')
            hold on
            xlim([min(x),max(x)])
            subplot(2,4,2)
            %figure('Name','Stoichiometry Distribution','NumberTitle','off')
            [plotH,~,x]=KDFplotH(trackArray(trackArray(:,8)>0,12),10);
            set(plotH,'Color',colourmap)
            xlabel('Spot Size Cell Parallel (nm)')
            ylabel('Probability Density')
            title('Linked Trajectories')
            hold on
            xlim([min(x),max(x)])
            subplot(2,4,3)
            %figure('Name','Diffusion Const Distribution','NumberTitle','off')
            [plotH,~,~]=KDFplotH(trackArray(trackArray(:,8)==0,13),10);
            set(plotH,'Color',colourmap)
            xlabel('Spot Size Cell Perpendicular (nm)')
            ylabel('Probability Density')
            title('Unlinked Trajectories')
            %xlim([-0.2,0.6])
            hold on
            subplot(2,4,4)
            %figure('Name','Diffusion Const Distribution','NumberTitle','off')
            [plotH,~,~]=KDFplotH(trackArray(trackArray(:,8)>0,13),10);
            set(plotH,'Color',colourmap)
            xlabel('Spot Size Cell Perpendicular (nm)')
            ylabel('Probability Density')
            title('linked Trajectories')
            %xlim([-0.2,0.6])
            hold on
            subplot(2,4,5)
            title('unlinked Trajectories')
            scatter(trackArray(trackArray(:,8)==0,2)/Isingle,trackArray(trackArray(:,8)==0,12),50,colourmap,'filled');
            hold on
            ylabel('Spot Size Cell Parallel (nm)')
            xlabel('Stoichiometry')
            hold on
            subplot(2,4,6)
            title('linked Trajectories')
            scatter(trackArray(trackArray(:,8)>0,2)/Isingle,trackArray(trackArray(:,8)>0,12),50,colourmap,'filled');
            hold on
            ylabel('Spot Size Cell Parallel (nm)')
            xlabel('Stoichiometry')
            subplot(2,4,7)
            title('unlinked Trajectories')
            scatter(trackArray(trackArray(:,8)==0,2)/Isingle,trackArray(trackArray(:,8)==0,13),50,colourmap,'filled');
            hold on
            ylabel('Spot Size Cell Perpendicular (nm)')
            xlabel('Stoichiometry')
            hold on
            subplot(2,4,8)
            title('linked Trajectories')
            scatter(trackArray(trackArray(:,8)>0,2)/Isingle,trackArray(trackArray(:,8)>0,13),50,colourmap,'filled');
            hold on
            ylabel('Spot Size Cell Perpendicular (nm)')
            xlabel('Stoichiometry')
        case 9
            subplot(2,3,1)
            scatter(stoichs,trackArray(:,12),'.')
            hold on
            AvData=binData(stoichs,trackArray(:,12),2,0:2:max(stoichs));
            hold on
            title('All tracks')
            ylabel('Spot Size Cell Parallel (nm)')
            xlabel('Stoichiometry')
            subplot(2,3,2)
            scatter(trackArray(trackArray(:,8)>0,2)/Isingle,trackArray(trackArray(:,8)>0,12),'.')
            hold on
            AvData=binData(trackArray(trackArray(:,8)>0,2)/Isingle,trackArray(trackArray(:,8)>0,12),2,0:2:max(stoichs));
            hold on
            title('Linked tracks')
            ylabel('Spot Size Cell Parallel (nm)')
            xlabel('Stoichiometry')
            subplot(2,3,3)
            scatter(trackArray(trackArray(:,8)==0,2)/Isingle,trackArray(trackArray(:,8)==0,12),'.')
            hold on
            AvData=binData(trackArray(trackArray(:,8)==0,2)/Isingle,trackArray(trackArray(:,8)==0,12),2,0:2:max(stoichs));
            hold on
            title('Unlinked tracks')
            ylabel('Spot Size Cell Parallel (nm)')
            xlabel('Stoichiometry')
            
            
            subplot(2,3,4)
            KDFplot(trackArray(:,12)./(stoichs));
            
            title('All tracks')
            ylabel('probability')
            xlabel('spacing/molecule (nm)')
            subplot(2,3,5)
            KDFplot(trackArray(trackArray(:,8)>0,12)./(trackArray(trackArray(:,8)>0,2)/Isingle));
            
            
            title('Linked tracks')
            ylabel('probability')
            xlabel('spacing/molecule (nm)')
            subplot(2,3,6)
            KDFplot(trackArray(trackArray(:,8)==0,12)./(trackArray(trackArray(:,8)==0,2)/Isingle));
            
            
            title('Unlinked tracks')
            ylabel('probability')
            xlabel('spacing/molecule (nm)')
            
        case 10
            colocClust=[];
            clustMod=[];
            coincidentMod=[];
            for c=1:max(AllSpots2linked(:,17)) % loop over all cells analysed
                possSpots=find(AllSpots2linked(:,17)==c);
                for s=min(possSpots):max(possSpots)
                    %    try
                    colocClustTemp=[];
                    colocClustTemp(:,1)=((AllSpots1linked(AllSpots1linked(:,17)==c,1)-AllSpots2linked(s,1)).^2.+...
                        (AllSpots1linked(AllSpots1linked(:,17)==c,2)-AllSpots2linked(s,2)).^2).^0.5;
                    colocClustTemp(:,2)=AllSpots1linked(AllSpots1linked(:,17)==c,5)/params.IsingleCh1;
                    colocClustTemp(:,3)=c;
                    
                    [row,col]=find(segArray(:,:,c)); %row with spots(:,2)
                    randInd=randsample(length(row),1);
                    %  randCoord(1)=range(AllSpots1linked(AllSpots1linked(:,17)==c,1))*rand(1)+min(AllSpots1linked(AllSpots1linked(:,17)==c,1))
                    colocClustTemp(:,4)=((AllSpots1linked(AllSpots1linked(:,17)==c,1)-col(randInd)).^2.+...
                        (AllSpots1linked(AllSpots1linked(:,17)==c,2)-row(randInd)).^2).^0.5;
                    colocClust=cat(1,colocClust,colocClustTemp);
                    %                     catch
                    %                     end
                    for r=1:length(row)
                        clustModtemp=((AllSpots1linked(AllSpots1linked(:,17)==c,1)-col(r)).^2.+...
                            (AllSpots1linked(AllSpots1linked(:,17)==c,2)-row(r)).^2).^0.5;
                        %  coincidentModTemp=rand(1)*params.colocError;
                        coincidentModTemp=normrnd(params.colocError/2,params.colocError/2);
                        clustMod=cat(1,clustMod,clustModtemp);
                        coincidentMod=cat(1,coincidentMod,coincidentModTemp);
                    end
                end
                
                
            end
            %             length(colocClust),
            %             length(clustMod),
            [counts,x]=histcounts(colocClust(:,1),0:0.5:max(colocClust(:,1)));
            PwD=pdist(AllSpots1linked(:,1:2));
            [counts3,x3]=histcounts(PwD,0:0.5:max(PwD));
            [counts2,x2]=histcounts(clustMod,0:0.5:max(clustMod));
            [counts4,x4]=histcounts(coincidentMod,0:0.5:max(clustMod));
            counts3=(counts3/sum(counts3))*sum(counts);
            counts2=(counts2/sum(counts2))*sum(counts);
            %  x4(counts4==max(counts4))
            counts4=(counts4/max(counts4))*counts(counts4==max(counts4));
            %  subplot(1,2,1)
            %   bar((0.25:0.5:(max(PwD)-0.25))*params.pixelSize,counts2)
            %hold on
            ylabel('Distance from Green spots')
            ylim([0,max([counts,counts2])])
            
            subplot(1,3,1)
            bar((0.25:0.5:(max(colocClust(:,1))-0.25))*params.pixelSize,counts)
            hold on
            %           plot((0.25:0.5:(max(colocClust(:,4))-0.25))*params.pixelSize,counts2)
            plot((0.25:0.5:(max(PwD)-0.25))*params.pixelSize,counts3,'LineWidth',2)
            plot((0.25:0.5:(max(clustMod)-0.25))*params.pixelSize,counts2,'LineWidth',2)
            xlabel('Distance from Red spots')
            ylim([0,max([counts,counts2])])
            legend('Red spots','Green pairwise distance','Random')
            subplot(1,3,2)
            bar((0.25:0.5:(max(colocClust(:,1))-0.25))*params.pixelSize,counts)
            hold on
            %           plot((0.25:0.5:(max(colocClust(:,4))-0.25))*params.pixelSize,counts2)
            plot((0.25:0.5:(max(PwD)-0.25))*params.pixelSize,counts3,'LineWidth',2)
            plot((0.25:0.5:(max(clustMod)-0.25))*params.pixelSize,counts2,'LineWidth',2)
            plot((0.25:0.5:(max(clustMod)-0.25))*params.pixelSize,counts4,'LineWidth',2)
            xlabel('Distance from Red spots (microns)')
            ylim([0,max([counts,counts2])])
            xlim([0,1])
            legend('Red spots','Green pairwise distance','Random','Colocalised with red model')
            subplot(1,3,3)
            plot((0.25:0.5:(max(colocClust(:,1))-0.25))*params.pixelSize,cumsum(counts/sum(counts)),'LineWidth',2)
            hold on
            plot((0.25:0.5:(max(PwD)-0.25))*params.pixelSize,cumsum(counts3/sum(counts3)),'LineWidth',2)
            plot((0.25:0.5:(max(clustMod)-0.25))*params.pixelSize,cumsum(counts2/sum(counts2)),'LineWidth',2)
            
            xlabel('Distance from Red spots')
            legend('Red spots','Green pairwise distance','Random')
            
            sum(counts(x<5))/sum(counts3(x3<5)),
            break
            
        otherwise
    end
    
    cd(DataDir)
    
    %%SAVE ANALYSED RESULTS in MAT structure
    if params.saveOutput==1
        cd ..
        cd ..
        %analysisname = 'category'
        %analysisfile = strcat(params.DataDir,'\',analysisname,'_output.mat'); %with autogenerated name
        analysisfile = strcat('output_plot',num2str(params.plotSelect),'.mat'); %simple name
        save(analysisfile,'output','params');
        disp('Copy info saved.')
        cd(DataDir)
    end
end
