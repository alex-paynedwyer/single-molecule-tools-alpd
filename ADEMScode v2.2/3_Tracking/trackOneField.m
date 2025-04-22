function [SpotsCh1, SpotsCh2, frame_average,p, metadata, image_data,spotImages] = trackOneField(fileName,p)

% Code for tracking bright foci in image stacks
%INPUTS
%filename: name of image stack or folder containing series of tifs
%p: parameter structure (see below for details) if not set default values
%used
%OUTPUTs
% SpotsCh1/2 is an array, each row contains the information for a foci found in an image frame in the series. The columns contain the following information:
% 
% 1. X coordinate (pixels)
% 2. Y coordinate (pixels)
% 3. Clipping_flag (a switch, please ignore)
% 4. Mean local background pixel intensity
% 5. Total foci intensity, background corrected
% 6. X foci sigma width
% 7. Y foci sigma width
% 8. foci central intensity (peak intensity in a fitted Gaussian)
% 9. Frame number the foci was found in
% 10. Trajectory number, foci in the same trajectory have the same trajectory number
% 11. Signal to noise ratio
% 12. Frame in which laser exposure began
% 
%frame_average: frame average image
%p: output parameters
% meta_data: any meta_data stored in the image
%image_data: the raw imaging data tracked
% spotImages: every foci tracked - NOT RECOMMENDED TO USE

%%
% NEW IN 2.9 (ALPD)
% Opens CZI directly and performs maximum intensity projection in Z for 2D tracking of confocal volumes
% Checks for saturated frames
% Able to crop out specified left and right margins of the field
% 
% NEW IN 2.81
% Rsquare of 2D gaussian fits;
% RSQUARE REPLACES CLIPPING FLAG IN foci ARRAY
% NEW IN 2.5
% can fit a fully rotating Gaussian using GaussSwitch=3;
% THETA REPLACES CLIPPING FLAG IN foci ARRAY
% NEW IN 2.4
% - choose to run in parallel or not by setting p.useParallel=1;
%- >2x speed improvements from optomisation in linking trajectories and
%removing coincident ones
%-2 levels of output set with show_output and show_all_output
% NEW IN 2.0
% -set useBioFormats to open ANY IMAGE FORMAT and extract metadata
% NB need file extension for this mode, NB EXCEPT SIF FILES?
% -removes candidates which are too close together to produce distinct
% spots, saves time
% -new iterative Gaussian masking to determine PSFwidth but only 1D, gives
% ~8x speed increase, set GaussSwitch=1 for this or =2 for full 2D fitting
% NEW IN 1.3
% ExtractImageSequence3 can now read in ascii data
% If FramesToTrack==0, then code tracks all frames available

%%

% Final code for opening tif data, calculating frame average, using this to determine
%cell boundary, then looping over user defined frames, finding spots,
%identifying their centres and accepting them if they meet criteria, then
%linking these foci together into trajectories

%readData=1 if reading tif file or 0 if already loaded

if exist('p','var')==1
    %Read in parameters
else
    % Number of frames to average over; used to choose foci in cursor mode and output
    p.noFrames=5; % number of frames to average over, starting from start_frame. (default: 5, end_frame-start_frame)
    
	% PARAMETERS for finding foci centres
    % The total integrated foci intensity is a bgnd corrected one, inside a circular mask of radius inner_circle_radius.
    p.subarray_halfwidth = 8; % (Default: 8 pixels). Halfwidth of image square subarray
    % ROI, typically a square of size 17x17 pixels.
    p.inner_circle_radius = 5; % (Default: 5 pixels). Radius of inner circular mask that moves inside the fixed square subarray.
    p.gauss_mask_sigma = 2; % (Default: 2 pixels). Size in pixels of the applied Gaussian mask.
    p.guess_sigma_Fit = 3; % starting guess for Gaussian fit of brightfoci intensity (default = 3).
    
    % PARAMETERS for deciding if we accept a foci centre
    p.sigmaFit_min = 0;  % minimum acceptable sigma of gaussian fit to spot, in pixels (2) (-3).
    p.sigmaFit_max = p.inner_circle_radius; % maximum acceptable sigma of gaussian fit to spot, in pixels (4) (3).
    p.SNR_min = 0.3;%0.4; % minimum acceptable signal-to-noise ratio (Usually set to level at which foci are ~50% false positives, SNR 0.3 for Prime95B or SNR 0.4 for Andor iXon camera)
    
    % PARAMETERS for eliminating coincident spots:
    p.d_min=1; % distance (in pixels) for eliminating coincidences
    
    % PARAMETERS for building trajectories:
    % For linking foci in current and previous frames:
    p.d_01_max = 8; % max distance in pixels between foci centres in current and previous frames, for linking them into a trajectory (8)
    p.Iratio_01_min = 0.5; % min ratio of total foci intensities (after bgnd subtraction) (0.5).
    p.Iratio_01_max = 3; % max ratio of total foci intensities (after bgnd subtraction) (frame k-1/frame k) (large enough value (3) to account for blinking).
    p.SigmaRatio_01_min = 0.5; % min ratio of foci widths (sigma of Gaussian fit) (0.5).
    p.SigmaRatio_01_max = 3; % max ratio of foci width (sigma of Gaussian fit) (2).
    p.error_set=0.05; %error in iterative Gaussian masking
    p.exclude_region = 0; % Parameter to exclude a slit/overlap region from the middle if Csplit>0;
    p.exclude_edge_left = 0; % Parameter to exclude slit regions from the sides if Csplit>0;
    p.exclude_edge_right = 0; % Parameter to exclude slit regions from the sides if Csplit>0;
    p.start_channel=1; %Channels to use
    p.end_channel=2;
    p.disk_radius = 5; % for finding foci in image
    p.ALEX=0; % Switch, if ALEX experiment=1
        	
    % There are 2 methods for finding candidate foci which work better for
    % different datasets or if set =3 runs both and uses all spots
    % (recommended)
    p.CandidateFindMethod=3;
    
    % Choose to iterate a Gaussian to determine PSF width (1D only) ==1
    % or explicitly fit a 2D Gaussian using MATLAB fitting routines==2
    % or fit a fully rotating Gaussian==3
    p.GaussSwitch=1;
    % switch to use parallel processing or not
    p.useParallel=1;
    %Distance to remove candidates which are too close together, set to 0 to ignore
    p.Candidate_d_min=0;
    
    %set to one and saves an array with every foci image saved
    p.spotImageSave=0;
    %gaussian=1 if running a gaussian filter over image data before finding foci
    p.gaussian=0;
    
    % Use this to specify interesting foci with the cursor, rather than autodetecting them
    p.use_cursor=0;
    
    % Use this to open files with BioFormats to open non tiff files
    p.useBioFormats=0;
    p.CZI=0;
    p.print_metadata=0;
    p.bitDepth = 12;
    
    p.all=1; %Use this keyword to load entire image file
    p.startFrame=1; % Or specify start and end frames if p.all=0
    p.endFrame=300;
    %Specify how many frames to track after the laser has switched on
    p.FramesToTrack=0;
    
    % Set this to 1 if there are blank frames before laser turns on or shutter opens
    p.DetermineFirstFrames=0;
    
    % Switch, =1 to determine laser on time with differential rather than max intensity
    p.use_diff=0;
    % CSplit defines how the channels are split, =0 for whole frame (no split), 1 for left/right and 2 for up/down
    p.CSplit=1;
    %If this =1, then graphs will appear
    p.show_output=0;
    % set this and graphs will appear for every spot!
    p.show_all_output=0;
    p.show_text_output=1;
    
end
%Initialise foci variables
SpotsCh1=[];
SpotsCh2=[];

%% CREATE FITTYPE AND OPTIONS
switch p.GaussSwitch
    case 1
        % don't need any fitting stuff
        myfit=[];
        options=[];
    case 2
        % Create fit type for constrained 2D Gaussian fit to spots
        myfit = fittype('Ibg_avg+(Isp./(2.*pi.*sdx.*sdy))*exp(-(((x-x0).^2)./(2.*sdx^2)+((y-y0).^2)./(2.*sdy^2)))',...
            'problem', {'Ibg_avg','Isp','x0','y0'}, 'independent', {'x', 'y'}, 'dependent', 'z');
        % Fit options:
        options = fitoptions(myfit);
        options.StartPoint = [p.guess_sigma_Fit, p.guess_sigma_Fit];
        options.Lower = [p.sigmaFit_min, p.sigmaFit_min];
        options.upper = [p.sigmaFit_max, p.sigmaFit_max];
    case 3
        myfit = fittype('Ibg_avg+ (Isp./(2.*pi.*sdx.*sdy))*exp( - ((cos(theta)^2/(2*sdx^2) + sin(theta)^2/(2*sdy^2))*(x-x0).^2 - 2*(-sin(2*theta)/(4*sdx^2) + sin(2*theta)/(4*sdy^2))*(x-x0).*(y-y0) + (sin(theta)^2/(2*sdx^2) + cos(theta)^2/(2*sdy^2))*(y-y0).^2))',...
            'problem', {'Ibg_avg','Isp','x0','y0'}, 'independent', {'x', 'y'}, 'dependent', 'z');
        options = fitoptions(myfit);
        options.StartPoint = [p.guess_sigma_Fit, p.guess_sigma_Fit,0];
        options.Lower = [p.sigmaFit_min, p.sigmaFit_min,0];
        options.upper = [p.sigmaFit_max, p.sigmaFit_max,2*pi];
end


%% OPEN DATA

%Open tif data from image_label which is whatever string is before the file
%extension
if p.useBioFormats==1
    if isa(fileName,'char')==1
        if p.CZI==1
            if p.all==1
                [image_data,metadata,image_X,image_Y,numFrames,numSlices,voxelUnit,voxelSizeX,voxelSizeY,voxelSizeZ]=openCZI_v1_1(fileName,p.print_metadata);
            else
                [image_data,metadata,image_X,image_Y,numFrames,numSlices,voxelUnit,voxelSizeX,voxelSizeY,voxelSizeZ]=openCZI_v1_1(fileName,p.print_metadata,p.startFrame,p.endFrame);
            end
            image_data=max(image_data,[],4); %max Z projection
        else
            if p.all==1
                [image_data, metadata]=imEx1(fileName);
            else
                [image_data, metadata]=imEx1(fileName,p.startFrame,p.endFrame);
            end
            [image_Y,image_X,numFrames]=size(image_data);
        end
    else
        image_data=fileName;
        [image_Y,image_X,numFrames]=size(image_data);
    end
else
    if isa(fileName,'char')==1
        [numFrames, image_Y, image_X, image_data, ~] = ExtractImageSequence3(fileName, p.all, p.startFrame, p.endFrame);
        metadata='sorry metadata not available in this mode';
    else
        image_data=fileName;
        [image_Y,image_X,numFrames]=size(image_data);
        metadata='sorry metadata not available in this mode';
    end
end

disp(['Dimensions: ', num2str(image_X), ' x ', num2str(image_Y)])
disp(['NumFrames: ', num2str(numFrames)])
%[numFrames, frame_Ysize, frame_Xsize, image_data, image_path] = extract_image_sequence_dataAWarray(image_label, all);
disp('data loaded')

%% ROTATE CHANNELS IF HORIZONTALLY SPLIT
if p.CSplit==2
    disp('Horizontal camera split; rotating stack clockwise.')
    image_data = rot90(image_data,3);
end

%% DETERMINE LASER ON FRAME

%Determine when laser turned on, detects both channels separately if ALEX
if p.DetermineFirstFrames==1
    [firstLeft, firstRight, ~, ~] = LaserOn3(image_data, p.use_diff,p.ALEX);
    disp('start determined')
else
    if p.ALEX==1
        if p.leftfirst==1
            firstLeft=1; firstRight=2;
        else
            firstLeft=2; firstRight=1;
        end
        if p.darkfirst==1
            firstLeft=firstLeft+1; firstRight=firstRight+1;
        end        
    else
        firstLeft=1;
        firstRight=1;
    end
    
    if p.FramesToTrack==0
        p.FramesToTrack=numFrames-firstLeft;
    end
    
    %% CALCULATE FRAME AVERAGE
    
    %Calculate Frame Average of the data
    try
    if firstLeft<(size(image_data,3)-5)
        frame_average = FrameAverage3(image_data, p.noFrames, firstLeft,firstRight,p.ALEX);
        disp('average calculated')
    else
        disp('WARNING FRAME AVERAGE NOT CALCULATED AS STARTFRAME TOO CLOSE TO ENDFRAME')
        frame_average=image_data(:,:,firstLeft);
    end
    catch err
        disp('calculating frame average failed because:')
        disp(err.message)
    end
    
    %%
    for Ch=p.start_channel:p.end_channel
        % Initialise foci array
        spots=[];
        spotImages=[];
        
        
        if p.ALEX==0
            if Ch==1
                startFrame=firstLeft;
                endFrame=firstLeft+p.FramesToTrack;
            else
                startFrame=firstRight;
                endFrame=firstRight+p.FramesToTrack;
            end
            FrameInt=1;
        else
            FrameInt=2;
            if Ch==1
                startFrame=firstLeft;
                endFrame=firstLeft+p.FramesToTrack;
            else
                startFrame=firstRight;
                endFrame=firstRight+p.FramesToTrack;
            end
        end
        if endFrame>size(image_data,3)
            endFrame=size(image_data,3);
        end
        
        %Loop over frames
        for i=startFrame:FrameInt:endFrame
            if p.show_text_output==0
                %    h=waitbar(i-startFrame/p.FramesToTrack);
            end
            disp(strcat('Tracking frame: ',num2str(i)))
            %% FIND CANDIDATE SPOTS
            %% Divide image into two channels, left and right
            if Ch==1
                if p.show_text_output==1
                    disp('Ch1')
                end
                switch p.CSplit
                    case 0 %full frame
                        frame=image_data(:,:,i);
                    case 1 %left slit window
                        %frame=image_data(:,:,i)-image_data(:,:,i);
                        frame=image_data(:,1+p.exclude_edge_left:round(size(image_data,2)/2-p.exclude_region),i);
                    case 2 %lower slit window becomes left slit window after rotation
                        %frame=image_data(1:round(size(image_data,2)/2-p.exclude_region),:,i)
                        %different clause not needed provided data rotation has
                        %worked properly.
                        frame=image_data(:,1+p.exclude_edge_left:round(size(image_data,2)/2-p.exclude_region),i);
                end
                SpotsCh1=[];
            elseif Ch==2
                if p.show_text_output==1
                    disp('Ch2')
                end
                switch p.CSplit
                    case 0  %full frame
                        frame=image_data(:,:,i);
                    case 1  %right slit window
                        frame=image_data(:,round(size(image_data,2)/2+p.exclude_region)+1:end-p.exclude_edge_right,i);
                    case 2  %upper slit window becomes right slit window after rotation
                        %frame=image_data(round(size(image_data,2)/2+p.exclude_region):end,:,i);
                        %different clause not needed provided data rotation has
                        %worked properly
                        frame=image_data(:,round(size(image_data,2)/2+p.exclude_region)+1:end-p.exclude_edge_right,i);
                end
                
                SpotsCh2=[];
            end
            
            %% Check for saturation
            if max(max(frame))>(2^p.bitDepth)
                disp('Warning! Frame includes saturated pixels.')
                continue
            end            
            
            %% For cursor mode
            if p.use_cursor==1
                if i==startFrame
                    pause on
                    % Display frame average to choose spots
                    if Ch==1
                        if p.CSplit==0
                            frame_averageCH=frame_average;
                        elseif p.CSplit==1
                            frame_averageCH=frame_average(:,1:round(size(image_data,2)/(p.end_channel-p.start_channel+1)-p.exclude_region));
                        elseif p.CSplit==2
                            frame_averageCH=frame_average(:,1:round(size(image_data,2)/(p.end_channel-p.start_channel+1)-p.exclude_region));
                        end
                    elseif Ch==2
                        if p.CSplit==0
                            frame_averageCH=frame_average;
                        elseif p.CSplit==1
                            frame_averageCH=frame_average(:,round(size(image_data,2)/2+p.exclude_region):end);
                        elseif p.CSplit==2
                            frame_averageCH=frame_average(:,round(size(image_data,2)/2+p.exclude_region):end);
                        end
                    end
                    imshow(frame_averageCH,[],'InitialMagnification','fit')%HM modify magnification
                    title('click a foci and hold alt key to select multiple spots, push any key when finished')
                    datacursormode on
                    pause
                    dcm_obj = datacursormode(1);
                    info_struct = getCursorInfo(dcm_obj);
                    %Loop over foci chosen and pull out co-ordinates
                    for q=1:size(info_struct,2)
                        Spot_coords=info_struct(q).Position;
                        y_estimate(q,1)=Spot_coords(2);
                        x_estimate(q,1)=Spot_coords(1);
                    end
                    close all
                end
            else
                % Create matrix of 1s where foci might be
                % Now 3 methods for doing this
                switch p.CandidateFindMethod
                    case 1
                        [result] = findSpots2(frame,2,p.disk_radius,p.gaussian,0);
                    case 2
                        [result] = findSpots3(frame,2,p.disk_radius,p.gaussian,0);
                    case 3
                        [result1] = findSpots2(frame,2,p.disk_radius,p.gaussian,0);
                        [result2] = findSpots3(frame,2,p.disk_radius,p.gaussian,0);
                        result=result1+result2;
                        result(result>1)=1;
						resulthin=bwmorph(result,'thin',Inf);
                end
                % Convert those to foci co-ordinates
                [y_estimateTemp, x_estimateTemp]=ind2sub(size(resultthin), find(resultthin));
                %  [y_estimate, x_estimate]=ind2sub(size(result), find(result));
                % Remove candidates which are too close together to yield separate spots
                if p.Candidate_d_min>0
                    [y_estimate,x_estimate]=MergeCoincidentCandidates2(y_estimateTemp, x_estimateTemp, p.Candidate_d_min);
                else
                    y_estimate=y_estimateTemp;
                    x_estimate=x_estimateTemp;
                end
                
            end
            if p.show_text_output==1
                disp('candidates found')
            end
            %Plot the candidate spots
            if p.show_output==1
                imshow(frame, [],'InitialMagnification','fit')%HM modify magnification
                hold on
                plot(x_estimate, y_estimate, 'o')
                hold off
                title('candidate spots')
                pause
            end
            %% FIT TO foci AND REJECT
            %Loop over all found foci
            spots_temp=zeros(size(x_estimate,1),12);
            spotImageTemp=zeros(p.subarray_halfwidth*2+1,p.subarray_halfwidth*2+1,size(x_estimate,1));
            if p.useParallel==1 % use a parfor loop
                parfor j=1:size(x_estimate,1)
                    if p.use_cursor==1
                        trajNo=j;
                        clip_override=1;
                    else
                        trajNo=0;
                        clip_override=0;
                    end
                    % Iterative gaussian masking to determine foci centre
                    [x_centre, y_centre, clipping_flag, Ibg_avg, Isp, Idata, bg_noise_std, mask_pixels,noConvergenceFlag]= ...
                        findSpotCentre4(frame,x_estimate(j),y_estimate(j),p,clip_override);
                                       
                    % Calculate foci signal to noise ratio
                    snr1=Isp/(bg_noise_std*mask_pixels);
                    
                    % check if already foci with those co-ords here
                    % isSpotAlready=sum(((spots_temp(:,1)-x_centre).^2+(spots_temp(:,2)-y_centre).^2).^0.5<p.d_min);
                    isSpotAlready=0;
                    % If in cursor mode OR a foci was found and it didn't clip
                    if p.use_cursor==1 || noConvergenceFlag==0 && clipping_flag==0
                        % Only store foci with good enough snr
                        if p.use_cursor==1 || snr1>p.SNR_min && isSpotAlready<1
                            switch p.GaussSwitch
                                case 1
                                    [sdx, sdy, Icent] = iterate1DgaussianFixedCenter3(frame,Ibg_avg, Isp, x_centre, y_centre,p);
                                    spots_temp(j,:)=[x_centre, y_centre, clipping_flag, Ibg_avg, Isp, sdx, sdy, Icent, i,trajNo, snr1, startFrame];
                                case 2
                                    [sdx, sdy, Icent,Rsquare] = fit2DgaussianFixedCenter4(frame,Ibg_avg, Isp, x_centre, y_centre,myfit,options,p);
                                    spots_temp(j,:)=[x_centre, y_centre, Rsquare, Ibg_avg, Isp, sdx, sdy, Icent, i,trajNo, snr1, startFrame];
                                case 3
                                    [sdx, sdy, Icent,theta] = fit2DRotGaussianFixedCenter3(frame,Ibg_avg, Isp, x_centre, y_centre,myfit,options,p);
                                    spots_temp(j,:)=[x_centre, y_centre, theta, Ibg_avg, Isp, sdx, sdy, Icent, i,trajNo, snr1, startFrame];
                                otherwise
                            end
                            % The foci array, 10th field is trajectory number,
                            % initialised to 0
                            spotImageTemp(:,:,j)=Idata;
                            
                        end
                    end
                    
                end
                spotImageTemp(:,:,spots_temp(:,1)==0)=[];
                spots_temp(spots_temp(:,1)==0,:)=[];
                % get rid of identical spots, within ~d_min of each other
                
                if isempty(spots_temp)==0
                    [~,ia,~]=uniquetol(spots_temp(:,1)+spots_temp(:,2),p.d_min/max(spots_temp(:,1)+spots_temp(:,2)));
                    %  [~,ia,~]=unique(spots_temp(:,1)+spots_temp(:,2));
                    %       try
                    spots_temp2=spots_temp(ia,:);
                    %               catch
                    %                   spots_temp(ia,:)
                    %                   spots_temp
                    %               end
                    spotImageTemp2=spotImageTemp(:,:,ia);
                end
                if exist('spots_temp2')
                    spots=cat(1,spots,spots_temp2);
                    spotImages=cat(3,spotImages,spotImageTemp2);
                end
            else %use a regular for loop
                for j=1:size(x_estimate,1)
                    if p.use_cursor==1
                        trajNo=j;
                        clip_override=1;
                    else
                        trajNo=0;
                        clip_override=0;
                    end
                    % Iterative gaussian masking to determine foci centre
                    [x_centre, y_centre, clipping_flag, Ibg_avg, Isp, Idata, bg_noise_std, mask_pixels,noConvergenceFlag]= ...
                        findSpotCentre3(frame,x_estimate(j),y_estimate(j),p,clip_override);
                    % Calculate foci signal to noise ratio
                    snr1=Isp/(bg_noise_std*mask_pixels);
                    
                    % check if already foci with those co-ords here
                    % isSpotAlready=sum(((spots_temp(:,1)-x_centre).^2+(spots_temp(:,2)-y_centre).^2).^0.5<p.d_min);
                    isSpotAlready=0;
                    % If in cursor mode OR a foci was found and it didn't clip
                    if p.use_cursor==1 || noConvergenceFlag==0 && clipping_flag==0
                        % Only store foci with good enough snr
                        if p.use_cursor==1 || snr1>p.SNR_min && isSpotAlready<1
                            switch p.GaussSwitch
                                case 1
                                    [sdx, sdy, Icent] = iterate1DgaussianFixedCenter3(frame,Ibg_avg, Isp, x_centre, y_centre,p);
                                    spots_temp(j,:)=[x_centre, y_centre, clipping_flag, Ibg_avg, Isp, sdx, sdy, Icent, i,trajNo, snr1, startFrame];
                                    
                                case 2
                                    [sdx, sdy, Icent,Rsquare] = fit2DgaussianFixedCenter4(frame,Ibg_avg, Isp, x_centre, y_centre,myfit,options,p);
                                    spots_temp(j,:)=[x_centre, y_centre, Rsquare, Ibg_avg, Isp, sdx, sdy, Icent, i,trajNo, snr1, startFrame];
                                    
                                case 3
                                    [sdx, sdy, Icent,theta] = fit2DRotGaussianFixedCenter3(frame,Ibg_avg, Isp, x_centre, y_centre,myfit,options,p);
                                    spots_temp(j,:)=[x_centre, y_centre, theta, Ibg_avg, Isp, sdx, sdy, Icent, i,trajNo, snr1, startFrame];
                                otherwise
                            end
                            % The foci array, 10th field is trajectory number,
                            % initialised to 0
                            
                            spotImageTemp(:,:,j)=Idata;
                            %      spots_temp(j,:)=[x_centre, y_centre, clipping_flag, Ibg_avg, Isp, sdx, sdy, Icent, i,trajNo, snr1, firstLeft];
                        end
                    end
                    
                end
                spotImageTemp(:,:,spots_temp(:,1)==0)=[];
                
                spots_temp(spots_temp(:,1)==0,:)=[]; %HM, gets rid of empty rows
                % get rid of identical spots, within ~d_min of each other
                
                if isempty(spots_temp)==0
                    [~,ia,~]=uniquetol(spots_temp(:,1)+spots_temp(:,2),p.d_min/max(spots_temp(:,1)+spots_temp(:,2)));
                    spots_temp2=spots_temp(ia,:);
                    spotImageTemp2=spotImageTemp(:,:,ia);
                end
                if exist('spots_temp2')
                    spots=cat(1,spots,spots_temp2);
                    spotImages=cat(3,spotImages,spotImageTemp2);
                end
            end
            
            
            if isempty(spots)==0
                if p.show_text_output==1
                    disp('Centres found and fitted')
                end
                
                %% PLOT ALL FOUND foci ON IMAGE
                %Plot all the found foci on image
                
                if p.show_output==1
                    imshow(frame,[],'InitialMagnification','fit') %HM modify magnification
                    hold on
                    title('Found ellipses on original image') %Plot ellipses on original image
                    if p.GaussSwitch==3
                        h=ellipse(spots(spots(:,9)==i,6),spots(spots(:,9)==i,7),spots(spots(:,9)==i,3),spots(spots(:,9)==i,1),spots(spots(:,9)==i,2),'b');
                    else
                        for k=min(find(spots(:,9)==i)):max(find(spots(:,9)==i))
                            text(spots(k,1)+3,spots(k,2)+3,num2str(k),'color','b')
                            rectangle('Position',[spots(k,1)-spots(k,6),spots(k,2)-spots(k,7),spots(k,6)*2,spots(k,7)*2],'Curvature',[1,1],'EdgeColor','b')
                        end
                    end
                    
                    hold off
                    pause
                end
                
                %% LINK foci INTO TRAJECTORIES
                if p.use_cursor==0
                    %Start at 2nd frame
                    if i>startFrame
                        %  if max(spots(:,9))==i
                        [spots]=LinkSpots4(spots, i, i-FrameInt, p);
                        %   end
                    end
                end
                if p.show_text_output==1
                    disp('trajectories determined')
                end
            else
                if p.show_text_output==1
                    disp('No foci found')
                end
            end
        end
        %% FINAL PLOT
        if isempty(spots)==0
            if p.show_output==1
                if max(spots(:,10))>0
                    subplot(2,3,2)
                    imshow(image_data(:,:,startFrame),[],'InitialMagnification','fit')%HM modify magnification
                    title('First frame with isolated foci superimposed')
                    hold on
                    plot(spots(spots(:,10)==0,1),spots(spots(:,10)==0,2),'o','color','r')
                    subplot(2,3,1)
                    imshow(image_data(:,:,startFrame),[],'InitialMagnification','fit')%HM modify magnification
                    title('First frame with tracks superimposed')
                    hold on
                    for i=1:max(spots(:,10))
                        traj_color=rand(1,3);
                        subplot(2,3,1)
                        title('Trajectory intensities vs frame num')
                        if Ch==1
                            plot(spots(spots(:,10)==i,1)+p.exclude_edge_left,spots(spots(:,10)==i,2),'-o','color',traj_color)
                        elseif Ch==2
                            plot(spots(spots(:,10)==i,1)+round(size(image_data,2)/2+p.exclude_region),spots(spots(:,10)==i,2),'-o','color',traj_color)
                        end
                        subplot(2,3,3)
                        hold on
                        plot(spots(spots(:,10)==i,9),spots(spots(:,10)==i,5),'-o','color',traj_color)
                        
                    end
                    title('Trajectory intensities vs frame num')
                    subplot(2,3,4)
                    
                    %    hist(spot_means)
                    hist(spots(:,5))
                    title('Histogram Trajectory Intensity Mean')
                    subplot(2,3,5)
                    
                    %                 hist(spot_sd)
                    %                 title('Histogram Trajectory Intensity SD')
                else
                    if p.show_text_output==1
                        disp('No trajectories found')
                    end
                end
            end
            %Assign foci to final variables
            if Ch==1
                if isempty(spots)==0
                    SpotsCh1=spots;
                    % Transform Channel1 so foci are in the right place on the uncropped image
                    SpotsCh1(:,1)=SpotsCh1(:,1)+p.exclude_edge_left;
                else
                    SpotsCh1=0;
                end
            elseif Ch==2
                if isempty(spots)==0
                    SpotsCh2=spots;
                    % Transform Channel2 so foci are in the right place on the uncropped image
                    SpotsCh2(:,1)=SpotsCh2(:,1)+round(size(image_data,2)/2+p.exclude_region);
                else
                    SpotsCh2=0;
                end
            end
            %  clear spots
        end
    end
    if p.show_text_output==0
        close(h)
    end
    
end