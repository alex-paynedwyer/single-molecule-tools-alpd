function [spotsBaselineCh1,spotsBaselineCh2]=overTracker(TifFile,DataFile,plotovertracks,selectinput,TifDir,DataDir)
noBaseFrames=10;
plotovertracks=0;
if selectinput==1
[TifFile,TifDir]=uigetfile('*.tif','Select the imaging data');
%cd(TifDir)
[DataFile,DataDir]=uigetfile('*.mat','Select the correct tracking data');
end
cd(DataDir);
load(DataFile);
subarray_halfwidth = p.subarray_halfwidth; % (Default: 8 pixels). Halfwidth of image square subarray
% ROI, typically a square of size 17x17 pixels.
inner_circle_radius = p.inner_circle_radius; % (Default: 5 pixels). Radius of inner circular mask that moves inside the fixed square subarray.
gauss_mask_sigma = p.gauss_mask_sigma; % (Default: 2 pixels). Size in pixels of the applied Gaussian mask.

cd(TifDir)
disp('loading image file')
[numFrames, frame_Ysize, frame_Xsize, image_data, image_path] = ExtractImageSequence3(TifFile(1:end-4), 1);

disp('initialising overtracking')
for ch=1:2
    switch ch
        case 1
            Spots=SpotsCh1;
            spotsBaselineCh1=Spots;
        case 2
            Spots=SpotsCh2; %need to use rawSpotsCh2 if tform applied!
            spotsBaselineCh2=Spots;
    end
    
    if ~isempty(Spots)
        %extract_image_sequence_dataAWarray(TifFile(1:end-4), 0, min(Spots(:,9)), max(Spots(:,9)));
        SpotLastFrame=zeros(max(Spots(:,10)),1);
        LFindex=zeros(max(Spots(:,10)),1);
        y_estimate=zeros(max(Spots(:,10)),1);
        x_estimate=zeros(max(Spots(:,10)),1);
        spot_num=size(Spots,1)+1;
        %Loop over trajectory numbers
        for k=1:max(Spots(:,10))
            [SpotLastFrame(k), LFindex(k)]=max(Spots(Spots(:,10)==k,9));
            SpotsX=Spots(Spots(:,10)==k,1);
            SpotsY=Spots(Spots(:,10)==k,2);
            y_estimate(k,1)=SpotsY(end);
            x_estimate(k,1)=SpotsX(end);
            %Loop over frames
            if SpotLastFrame(k)+noBaseFrames < max(Spots(:,9))
                for index=SpotLastFrame(k)+1:SpotLastFrame(k)+1+noBaseFrames
                    frame=image_data(:,:,index);
                    [x_centre, y_centre, clipping_flag, Ibg_avg, Isp, Idata, bg_noise_std, mask_pixels,noConvergenceFlag]= ...
                        findSpotCentre2noloop(frame,x_estimate(k),y_estimate(k),subarray_halfwidth,inner_circle_radius,gauss_mask_sigma,0.05, 0,0);
                    snr1=Isp/(bg_noise_std*mask_pixels);
                    switch ch
                        case 1
                            spotsBaselineCh1(spot_num,:)=[x_centre, y_centre, clipping_flag, Ibg_avg, Isp, inner_circle_radius, inner_circle_radius, -1, index,k, snr1, Spots(1,12)];
                        case 2
                            spotsBaselineCh2(spot_num,:)=[x_centre, y_centre, clipping_flag, Ibg_avg, Isp, inner_circle_radius, inner_circle_radius, -1, index,k, snr1, Spots(1,12)];
                    end
                    spot_num=spot_num+1;
                end
            end
        end
    end
end
cd(DataDir);
datafilename=strcat(DataFile(1:end-11),'_BASELINE.mat')
save(datafilename,'spotsBaselineCh1','spotsBaselineCh2');
disp('BASELINE file saved.')
%Plots the overtracked trajectories' brightness over time
%if flag is 0 plot nothing, if 1 or 2 plot one channel, if 3 plot both

if plotovertracks==1 || plotovertracks==3
figure
for i=1:max(spotsBaselineCh1(:,10))
plot(spotsBaselineCh1(spotsBaselineCh1(:,10)==i,9),spotsBaselineCh1(spotsBaselineCh1(:,10)==i,5))
hold on
end
hold off
end

if plotovertracks==2 || plotovertracks==3
figure
    for i=1:max(spotsBaselineCh1(:,10))
plot(spotsBaselineCh1(spotsBaselineCh1(:,10)==i,9),spotsBaselineCh1(spotsBaselineCh1(:,10)==i,5))
hold on
end
hold off
end

end