%% SCRIPT TO generate bead TRANSFORMS between dual channels, to be run before tracking single folders of data
%Use the output transform as 'p.beadtform' in the trackOneField/trackAllFields code.
%Run on the bead input tif image sequence (typically in the 'test' folder created when aligning the microscope)

%% Select bead image (must be same day and same region of interest/field of view as the dataset for tracking)

showtransformoutput = 1;
paramsbead.useavgproj = 0; %flag, if 0 then use first frame of bead video only, if 1 then use mean projection over all time frames; if 2 then use max projection over every 10th frame.
paramsbead.CSplit=1; %select camera split on the bead image (default =1 =left/right, alternative 2=top/bottom)
paramsbead.bitDepth=12; %default
paramsbead.topbottomcut = 0; % default=0, else use for cropping the edges of the bead image to fit the field of view
paramsbead.leftrightcut = 0; % default=0, else use for cropping the edges of the bead image to fit the field of view

[paramsbead.beadImageFile,paramsbead.beadImageDir]=uigetfile('*.tif');
cd(paramsbead.beadImageDir)
[beadnumFrames, beadimage_Y, beadimage_X, bead_image, ~] = ExtractImageSequence3(paramsbead.beadImageFile(1:end-4), 1, 1, 1);

if paramsbead.useavgproj==0
    bead_frame = bead_image(:,:,1); % first frame only
elseif paramsbead.useavgproj==1
    bead_frame = mean(bead_image,3); % 3 refers to time dimension
elseif paramsbead.useavgproj==2
    bead_frame = max(bead_image(:,:,1:10:beadnumFrames),[],3); % 3 refers to time dimension
end


gaussweight=0;
histeqflag=0;
[beadsframech1,beadsframech2,SlitSegmentObject]=prepRegImages(bead_frame,paramsbead,gaussweight,histeqflag);

%beadsframech1 = imadjust(bead_frame(:,1:beadimage_Y/2),[]); beadsframech2 = imadjust(bead_frame(:,beadimage_Y/2+1:beadimage_Y),[]);  

[xOff, yOff, mag, angle,beadtform]=ManualAlignBeadsTransform(beadsframech1,beadsframech2,showtransformoutput);

cropfilenameby = 21; %characters for '_MMStack_Pos0.ome.tif' standard on Micromanager 1.4
beadtformfilename = strcat(paramsbead.beadImageFile(1:end-cropfilenameby),'_TFORM.mat');

save(beadtformfilename,'beadtform','paramsbead');
cd ..;cd ..; % change the present working directory up two levels back to the strain folders.


%% Loop through nested folders to find image files
%for i=1:noDataDirs
    %if p.batchdirselect == 0
    %    DataDir=DataDirs(i)
    %end
    %cd(DataDir)
    
    % Directory to save the analysis data, duplicates hierarchy of data
    %ResultsPath=strcat(DataDir,'_ANALYSIS');
    %AltResultsPath=strcat(pwd,'\TrackingFiles');
    %mkdir(AltResultsPath)
    %mkdir(ResultsPath)