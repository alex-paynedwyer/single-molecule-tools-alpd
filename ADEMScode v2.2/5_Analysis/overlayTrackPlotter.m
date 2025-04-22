%plot tracks as overlay on a representative image

%% Load output file

showch1 = 1;
showch2 = 0;
plotallfovs = 1; %plots every track in an output file
plotanfov = 1; %index of the particular fov of interest;

disp('Please load TRACKS.mat');%select an FoV to visualise [manually select the same fov]
[trackName,trackDir]=uigetfile('*_TRACKS.mat');
cd(trackDir)
load(trackName); 

%% Plot relevant tracks from all fields as overlay

disp('Please load analysis output.mat');
[outputName,outputDir]=uigetfile('*.mat');
cd(outputDir)
load(outputName);

%finalimage = [load an appropriate image, eg. registered TIFF from FIJI, as array]
finalimage = frame_average;

%spots not tracks
%imshow(MAX_50nM(:,1:600),[]); hold on;
%snrvis=[0.5,0.5]; hold on;scatter(SpotsCh1(SpotsCh1(:,11)>snrvis(1),1),SpotsCh1(SpotsCh1(:,11)>snrvis(1),2)); hold on; scatter(SpotsCh2(SpotsCh2(:,11)>snrvis(2),1)-600,SpotsCh2(SpotsCh2(:,11)>snrvis(2),2));

%tracks overlay
numtracks=10000;
tracklengthmin=10;

brighttracks1=output.trackArrayCh1;
brighttracks1(1:numtracks,15)=1:numtracks;%1:length(brighttracks1);
if plotallfovs == 0
    fovnames=unique(brighttracks1(:,7));
    brighttracks1=brighttracks1(brighttracks1(:,7)==fovnames(plotanfov),:);
    
    %do stuff to select the correct image for the selected fov [not currently implemented]
    
end
brighttracks1=brighttracks1(brighttracks1(:,2)>0,:);
brighttracks1=brighttracks1(brighttracks1(:,14)>tracklengthmin,:);
brighttracks1list=brighttracks1(:,15);
imshow(finalimage(:,:),[]); hold on
if showch1 ==1
for i=1:length(brighttracks1list)
trackno=i; line(SpotsCh1(SpotsCh1(:,10)==trackno,1),SpotsCh1(SpotsCh1(:,10)==trackno,2));
end
end
if showch2==1
brighttracks2=output.trackArrayCh2;
brighttracks2(:,15)=1:length(brighttracks2);
if plotallfovs == 0
    fovnames=unique(brighttracks2(:,7));
    brighttracks2=brighttracks2(brighttracks2(:,7)==fovnames(plotanfov),:);
    
    %do stuff to select the correct image for the selected fov [not currently implemented]
end
brighttracks2=brighttracks2(brighttracks2(:,2)>200,:);
brighttracks2=brighttracks2(brighttracks2(:,14)>3,:);
brighttracks2list=brighttracks2(:,15);
for i=1:length(brighttracks2list)
trackno=i; line(SpotsCh2(SpotsCh2(:,10)==trackno,1)-600,SpotsCh2(SpotsCh2(:,10)==trackno,2));
end
end   
  