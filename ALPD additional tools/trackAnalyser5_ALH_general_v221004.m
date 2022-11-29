
function [trackArray,spotsInTracks,segArray,cellNo]=trackAnalyser5_ALH_general(spots,CellObject,previousTrackArray,frameLimit,frameTime,pixelSize,tracksFile,PSFsize,previousSegArray,cellNo,lowertracklimit,uppertracklimit)
%% function trackAnalyser
% does what you'd think, pulls any numbers from string in tracksFile and
% puts in column 7 of the array
%if PSFsize=0, no correction for PSF size or motion blur
trackArray=previousTrackArray;
segArray=previousSegArray;
spotsInTracks=[];

%switch motion blur correction on/off
if PSFsize>0
    sizeSwitch=1000;
else
    sizeSwitch=0;
end

if isempty(previousTrackArray)
  %  cellNo=0;
    spotNo=0;
else
  %  cellNo=max(previousTrackArray(:,1));
    spotNo=length(previousTrackArray(:,6));
end

for c=1:size(CellObject,3)
    if sum(CellObject(:,:,c),'all')>0
    cellCoord=[];
    [cellCoord(:,2), cellCoord(:,1)]=find(CellObject(:,:,c));
    cellOrientation=regionprops(CellObject(:,:,c),'orientation');
    spotInd=ismember(round(spots(:,1:2)),cellCoord,'rows');
    trajNo=unique(spots(spotInd,10));
    trajNo(trajNo==0)=[];
    if isempty(trajNo)
        continue
    end
    cellNo=cellNo+1;
    segArray(:,:,cellNo)=CellObject(:,:,c);
    for trajInd=1:length(trajNo)
        t=trajNo(trajInd);
        if min(spots(spots(:,10)==t,9))<(frameLimit+min(spots(:,9)))
            lentrack=length(spots(spots(:,10)==t,9));
            if lentrack>=lowertracklimit           %lentrack>3          %Crucial filter for two and/or three-spot tracks; this version also allows three-spot tracks if tracklimit set to 3 instead of 4
                if lentrack == lowertracklimit
                    tracklim = lowertracklimit;
                elseif lentrack>=uppertracklimit
                    tracklim = uppertracklimit;
                else
                    tracklim = lentrack;
                end
                spotNo=spotNo+1;
                SpotX=spots(spots(:,10)==t,1);
                SpotY=spots(spots(:,10)==t,2);
                SpotThet=spots(spots(:,10)==t,3)+(cellOrientation.Orientation./360*2*pi);
                SpotSX=spots(spots(:,10)==t,7);
                SpotSY=spots(spots(:,10)==t,6);
                SpotT=spots(spots(:,10)==t,9)-min(spots(spots(:,10)==t,9));
                SpotI=spots(spots(:,10)==t,5);
               
     
                % determine stoichiometry by fitting a line and using intercept
                pfit=polyfit(SpotT(1:tracklim)-1,SpotI(1:tracklim),1);
                %                 scatter(spots(spots(:,10)==t,9)-unique(spots(spots(:,10)==t,12)),spots(spots(:,10)==t,5))
                %                 hold on
                %                 plot(spots(spots(:,10)==t,9)-unique(spots(spots(:,10)==t,12)),(spots(spots(:,10)==t,9)-unique(spots(spots(:,10)==t,12)))*pfit(1)+pfit(2))
                %                 pause
                trackArray(spotNo,1)=cellNo; %cell ID number
                if pfit(2)>0 && pfit(1)<0
                    trackArray(spotNo,2)=pfit(2); %trajectory intensity (divide by Isingle for stoichiometry)
                else
                    trackArray(spotNo,2)=mean(SpotI(1:tracklim));
                end
                %diffusion coefficient
                [trackArray(spotNo,3),MSD,tau,LocPrecision]=getDiffusion3(spots(spots(:,10)==t,:),frameTime,pixelSize,tracklim);
                % first displacement D(1)
                trackArray(spotNo,4)=((SpotX(2)-SpotX(1)).^2+(SpotY(2)-SpotY(1)).^2)*(pixelSize).^2/(4*frameTime);
                trackArray(spotNo,5)=(SpotSX(1)^2+SpotSY(1)^2)^0.5; %spot size
                trackArray(spotNo,6)=t; %trajectory number
                trackArray(spotNo,7)= str2num(tracksFile(regexp(tracksFile,'\d')));
                trackArray(spotNo,8:11)=0; % these will be assigned later with colocalisation info
                % spot width along cell
                trackArray(spotNo,12)=max([abs(SpotSY(1)*cos(SpotThet(1))),abs(SpotSX(1)*cos(90-SpotThet(1)))])*pixelSize*1000-PSFsize-sizeSwitch*(2*frameTime*abs(trackArray(spotNo,3)))^0.5;
                % spot width perpendicular to cell
                trackArray(spotNo,13)=max([abs(SpotSY(1)*sin(SpotThet(1))),abs(SpotSX(1)*sin(90-SpotThet(1)))])*pixelSize*1000-PSFsize-sizeSwitch*(2*frameTime*abs(trackArray(spotNo,3)))^0.5;
               % spots(spots(:,10)==t,12)=cellNo;
               trackArray(spotNo,14)=lentrack;  %total number of spots in the track from one channel (approx. dwell time/interframe time)
               
               %update spot array with new information about tracks
               spotsInCell=spots(spots(:,10)==t,:);
               spotsInCell(:,17)=cellNo;
                spotsInTracks=cat(1,spotsInTracks,spotsInCell);
            end
        end
    end
    end
end
end


