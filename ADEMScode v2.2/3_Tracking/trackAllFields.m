%% Handy visualisation functions
%scatter(SpotsCh1(:,1),SpotsCh1(:,2)); hold on; scatter(SpotsCh2(:,1)-framewidth+0,SpotsCh2(:,2)+0);
%scatter(SpotsCh1(:,1),SpotsCh1(:,2)); hold on; scatter(SpotsCh2(:,1)-0,SpotsCh2(:,2)+0);

%% SET THE TRACKING PARAMETERS

% PARALLELISE the tracking computation
p.useParallel=1; %takes ~10 min to initialise but then runs faster over multiple cores.
p.batchdirselect = 1;  %using hardcoded folders below instead of selecting images folder to track at runtime
if p.batchdirselect == 1
    p.DataPrefix = "D:\MATLAB-current\Data\Arabidopsis\SlimVar\FLC colocalisation\Processed"
    p.DataFolders = [...
    "\NV",...
    "\V2W",...
    "\V6W",...
    "\V6WT14",...
    ];
end

% PARAMETERS for screening candidate foci
p.noFrames=5; % number of frames to average over, starting from start_frame. (default: 5, end_frame-start_frame)
p.subarray_halfwidth = 8; % (Default: 8 pixels). Halfwidth of image square subarray
p.inner_circle_radius = 5; % (Default: 5 pixels). Radius of inner circular mask that moves inside the fixed square subarray.
p.gauss_mask_sigma = 2; % (Default: 2 pixels). Size in pixels of the applied Gaussian mask.
p.guess_sigma_Fit = 3; % starting guess for Gaussian fit of local maximum intensity (default = 3).

% PARAMETERS for deciding if we accept local maxima as foci
p.SNR_min = 0.5;%0.2; % minimum possible signal-to-noise ratio (sifting)
p.sigmaFit_min = 0;  % minimum acceptable sigma of gaussian fit to spot, in pixels
p.sigmaFit_max = p.inner_circle_radius; % maximum acceptable sigma of gaussian fit to local maxima, in pixels
p.d_min=1; % distance (in pixels) for eliminating coincidences

% Camera setup
p.ALEX=0; % Switch, if ALEX experiment=1
% CSplit defines how the channels are split, =0 for whole frame (no split), 1 for left/right and 2 for up/down
p.CSplit=0;  % GFP/RFP and YFP/RFP modes, typically CSplit=1, Ch1 is RED and Ch2 is GREEN/YELLOW
p.registerchannels=0; %correct Ch2-Ch1 shift in Ch2 spots: (0) ignore; (1) using manual iteration; (2) using preexisting 'tform' transform from bead cursor mode, or (3) brightfield affine optimiser with bead transform 'tform' as starting point.
p.excludefirsttiff = 0;  % exclude the first tiff in the folder e.g. if it's brightfield

% Choosing what to track
p.allfields=1;  %Use this flag to load all fields of view
p.startField=27; %or curtail the number of fields to a certain range
p.endField=30;   

p.start_channel=1;   %Range usually set as 1,2. Set 1,1 or 2,2 if only one channel.
p.end_channel=1;

p.all=1; %Use this keyword to load entire image sequence
p.startFrame=1; % Or specify start and end frames if p.all=0; Ignored if p.all=1.
p.endFrame=30;
p.FramesToTrack=0; %Specify how many frames to track after the laser has switched on
% Set this to 1 if there are blank frames before laser turns on or shutter opens

p.leftfirst=0; %is the left channel triggered first? Usually not.
p.darkfirst=0; %is there a single dark frame at the start before any laser triggering?

p.DetermineFirstFrames=0;   %use to find the 'laser on' frame, e.g. you opened a shutter manually instead of the detector triggering the excitation
p.use_diff=0;  % Switch, =1 to determine laser on time with differential rather than max intensity

% PARAMETERS for building tracks & for colocalising foci in current and previous frames:
p.d_01_max = 8; %max distance in pixels between foci centres in current and previous frames, for linking them into a trajectory.
p.Iratio_01_min = 0.5; % min ratio of total foci intensities after background subtraction.
p.Iratio_01_max = 3; % max ratio of total foci intensities after background subtraction (frame k-1/frame k); large enough value (3) to account for blinking.
p.SigmaRatio_01_min = 0.5; % min ratio of foci widths (sigma of Gaussian fit).
p.SigmaRatio_01_max = 3; % max ratio of foci Gaussian width (sigma of  fit).
p.error_set=0.05; %error in iterative Gaussian masking
p.exclude_region = 0; % Parameter to exclude a slit/overlap region from the middle if Csplit>0;
p.exclude_edge_left = 0; % Parameter to exclude slit regions from the sides if Csplit>0;
p.exclude_edge_right = 0; % Parameter to exclude slit regions from the sides if Csplit>0;
p.disk_radius = 5; % for finding foci in image
p.topbottomcut = 0;
p.leftrightcut = 0; % setting to zero avoids additional automatic cropping in prepRegImages.

% Foci finding settings
p.CandidateFindMethod=3;
p.GaussSwitch=1;
p.Candidate_d_min=0;  %distance to remove candidates which are too close together
p.spotImageSave=0; %saves an array with each of the foci 
p.gaussian=0; %gaussian=1 if running a gaussian filter over image data before finding foci

% Options for loading unusual files
p.useBioFormats=0;
p.CZI=0; %used for CZI format files e.g. confocal images
p.bitDepth=12;

%Options for output
p.use_cursor=0;
p.print_metadata=0;
p.show_output=0;%If this =1, then graphs will appear at each step!
p.show_all_output=0; % set this and graphs will appear for each of the foci!
p.show_text_output=1;

SpotsCh1=[]; %Initialise Spot array variables
SpotsCh2=[];

%% Begin tracking script runtime: 
% Select dual colour bead transform as starting point for individual FOV refinements
if p.start_channel<=2 && p.end_channel>=2 && p.registerchannels>0
    if ~exist('p.beadtform')
       try  
            p.beadtform;
            disp('Initial transform found for dual channel registration: tform already in memory.')
       catch
            [beadtformFile,beadtformDir]=uigetfile('*_TFORM*','Select file containing initial bead transform');
            currentDir = pwd;
            cd(beadtformDir)
            load(beadtformFile);
            cd(currentDir)
            p.beadtform = beadtform;
            p.paramsbead = paramsbead;
            disp('Initial transform found for dual channel registration: tform from imported file.')
       end
    else
        disp('Initial transform found for dual channel registration in parameters.')
    end
else
    disp('Single channel tracking: no channel registration performed.')
end

%%Choose STRAIN directories
%Expected folder hierarchy = STRAIN/DATE/SAMPLE/FIELD  in format "*/2YYY-MM-DD/sample*/field*/*.tif"

if p.batchdirselect==1
    %%Hardcode a list of direct paths in full
    %DataDirs=[];
    %%Or hardcode the subdirectories
    DataPrefix = p.DataPrefix;
    DataFolders = p.DataFolders;
    DataDirs =[];
    for j = 1:length(DataFolders)
        DataDirs=[DataDirs,strcat(DataPrefix,DataFolders(j))];
    end
    noDataDirs = length(DataDirs);
else
    DataDir=uigetdir('*.mat','Select the strain folder to analyse')
    noDataDirs = 1;
end

%% Loop through nested folders to find image files
for i=1:noDataDirs
    if p.batchdirselect == 1
        DataDir=DataDirs(i)
    end
    cd(DataDir)

    % Directory to save the analysis data, duplicates hierarchy of data
    ResultsPath=strcat(DataDir,'_ANALYSIS');
    AltResultsPath=strcat(pwd,'\TrackingFiles');
    %mkdir(AltResultsPath)
    mkdir(ResultsPath)

    prebleachonly=0;  %Use if there's an extra prebleach fluorescence acquisition (a third video per folder)

    %directories
    DateDir=dir('*2*');  %date string must always include year '2YYY'.
    DateNum=size(DateDir);

    % Loop over different date folders
    for j=1:DateNum(1) %1:DateNum(1) %change as needed to exclude unwanted data
        cd(strcat(DateDir(j).folder,'\',DateDir(j).name))
        SampleDir=dir('*am*'); %find all subfolders with names containing 'Sample'
        SampleDir=SampleDir([SampleDir.isdir]>0);
        SampleNum=size(SampleDir);
        DatePath=strcat(ResultsPath, '\', DateDir(j).name);
        mkdir(DatePath)
        % Loop over sample folders
        for k=1:SampleNum(1)
            cd(SampleDir(k).name)
            CellDir=dir('*el*'); %find all subfolders with names containing 'field' or 'cell'
            CellDir = CellDir(~contains({CellDir.name},'DS_Store'));
            CellDir=CellDir([CellDir.isdir]>0);
            [CellDirSorted,cellindexlist]=natsort({CellDir.name});
            CellNum=size(CellDir);
            SamplePath=strcat(DatePath, '\', SampleDir(k).name);
            mkdir(SamplePath)
            if p.allfields==1
                p.startField=1;
                p.endField=CellNum(1);
            end
            % Loop over field folders containing tifs
            for l=p.startField:p.endField    %natural order (1,2,..10) not alphanumeric (1,10..,2)
                strcat('Tracking field ',num2str(l))
                BFindex=[];
                cd(CellDirSorted{l})
                 
                CellPath=strcat(SamplePath, '\', CellDirSorted{l});
                mkdir(CellPath)
                TifFiles=dir('*.tif');
                
                % Check there are some tifs to analyse
                EnoughTifs=size(TifFiles);
                if EnoughTifs(1)>0
                    %  Check there is a full acquisition by size of file
                    if max([TifFiles(:).bytes]) > 3005
                        pwd
                        
                        FLUORname=TifFiles(find([TifFiles.bytes]== max([TifFiles(:).bytes]))).name;

                        % Open 1 BF image
                        BFindex=find([TifFiles.bytes]== min([TifFiles(:).bytes]));
                        BFInfo=imfinfo(TifFiles(BFindex(end)).name);
                        BFImage=imread(TifFiles(BFindex(end)).name);
                        BFImage=mean(BFImage,3);
                        framewidth = size(BFImage,2)/2;
                        BFFileName=strcat(CellPath, '\', 'BF.tif');
                        BFFileName2=strcat(AltResultsPath, '\',TifFiles(1).name, 'BF.tif');


                        %% TRACK DATA in earnest

                        % Choose which images to track in each folder based on index 'b' (usually timestamp order)
                        if p.excludefirsttiff == 0
                            bmin=1;                    %INCLUDE FIRST TIFF
                        elseif p.excludefirsttiff == 1
                            bmin=2;                 %TO EXCLUDE BRIGHTFIELD
                        end
                        for b=bmin:length(TifFiles)
                            %for b=3:length(TifFiles)  %TO EXCLUDE BRIGHTFIELD and PREBLEACH
                            %prebleachonly=1; for b=2:2  %FOR PREBLEACH STEP ONLY

                            fileToTrack=TifFiles(b).name

                            [SpotsCh1, SpotsCh2, frame_average,~,~,~] = trackOneField(fileToTrack(1:end-4),p);

							%POST-TRACKING CHANNEL REGISTRATION
                            rawSpotsCh2 = SpotsCh2; %save the original spots as tracked
                            if p.registerchannels~=0&&p.start_channel<=2 && p.end_channel>=2
                                % PREPARE IMAGE REGISTRATION based on shift in brightfield of each FOV
                                [regframech1,regframech2,framech1,framech2,CellObject]=prepRegImages(BFImage,p,1,1); %crop images to each channel and weight by the central laser area
    
                                if p.registerchannels==1
                                    % APPLY INTERCHANNEL REGISTRATION MANUALLY
                                    %try to correct the overlap iteratively using circshift to translate channel 2 onto channel 1
                                    Ch2toCh1shift_x = +3;
                                    Ch2toCh1shift_y = -8;
                                    shiftedbfch2=circshift(regframech2,[Ch2toCh1shift_y,-Ch2toCh1shift_x]);
                                    figure;imshow(cat(3,regframech1,shiftedbfch2,regframech1),[])
                                    
                                elseif p.registerchannels==2
                                    % APPLY THE AUTOMATIC AFFINE TRANSFORM DIRECTLY FROM EXISTING INPUT i.e. BEADTFORM
                                    SpotsCh2(:,1)=SpotsCh2(:,1)-framewidth;
                                    alignedSpotsCh2 = applySpotsTform(SpotsCh1,SpotsCh2,p.beadtform,framewidth,0,0);
                                    tform = p.beadtform;
                                    SpotsCh2 = alignedSpotsCh2;  %apply the transform using 'transformPointsForward' and overwrite
                                    SpotsCh2(:,1)=SpotsCh2(:,1)+framewidth;
                                    disp('Bead-derived registration transform applied to channel 2 spots.')
                                    
                                elseif p.registerchannels==3
                                    % REFINE AFFINE TRANSFORM USING BRIGHTFIELD SIMILARITY OPTIMISER BEFORE APPLYING
                                    disp('Optimising bead-derived channel registration with brightfield images.')
                                    invertflag = 0;
                                    SpotsCh2(:,1)=SpotsCh2(:,1)-framewidth;
                                    [imageCh2start,imageCh2adj,tform]=transformRefinement(regframech1,regframech2,p.beadtform,invertflag,0);
                                    alignedSpotsCh2 = applySpotsTform(SpotsCh1,SpotsCh2,tform,framewidth,0,0);
                                    SpotsCh2 = alignedSpotsCh2;  %apply the transform using 'transformPointsForward' and overwrite
                                    SpotsCh2(:,1)=SpotsCh2(:,1)+framewidth;
                                    disp('Brightfield-optimised registration transform applied to channel 2 spots.')
                                end
                            else
                                tform=affine2d([1,0,0;0,1,0;0,0,1]);  %otherwise use blank, identity transform between channels
                                disp('One channel only - defaulted to blank transform.')
                            end
                            

                            % Save SpotsCh1 to spots.mat and frame average;
                            datafilename=strcat(CellPath, '\',fileToTrack(1:end-4),'_TRACKS.mat');
                            %datafilename=strcat(CellPath, '\',fileToTrack(1:end-16),'_TRACKS.mat'); %removes '_MMStack.ome' from analysed filenames
                            %altdatafilename=strcat(AltResultsPath,'\',fileToTrack(1:end-16),'_TRACKS.mat');


                            %% SAVE TRACKS.mat to ANALYSIS path

                            save(datafilename,'SpotsCh1','SpotsCh2','p','frame_average','tform');
                           
                            %segdatafilename=strcat(CellPath,'\',FLUORname(1:end-4),'_segmentation.mat');
                            %altsegdatafilename=strcat(AltResultsPath,'\',FLUORname(1:end-4),'_segmentation.mat');
                            %save(segdatafilename,'CellObject','Cellframe0','tform');
                            %save(altsegdatafilename,'CellObject','Cellframe0','tform');
                        end
                    end
                end
                cd ..
            end
            cd ..
        end
        cd ..
    end
end
%end
