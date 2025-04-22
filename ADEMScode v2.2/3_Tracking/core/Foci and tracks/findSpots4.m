function [result] = findSpots4(frame,disk_radius, output)
%
% Find candidate (fluorescence) bright spots on images (of cells).
% The candidate spots can later be used as input of function
% "findSpotCentre1frame".
%
%
% INPUTS:
% frame: Input image, a single frame from a recorded sequence.
%
% method: two different methods, either 1 or 2.
%
% gaussian: a switch which runs a guassian filter if set to 1
%
% disk radius in pixels for top hat transformation, quite critical to the method for finding candidate bright spots.
%
%
% output=1 to see figures and advance manually
%
%
%
% OUTPUTS:
%result is matrix same size as the image with 1s at candidate spots
%

if output==1
    pause on
end

frame0 = mat2gray(frame);
if output==1
    imshow(frame0,[])
    title('original image frame')
    pause
end


se = strel('disk',disk_radius); % structural element, disk of radius 5 pixels (default). The radius is quite critical to the method.
Signal = imtophat(frame0,se); % enhanced signal-image via "top hat" transformation (evens out the background).
if output==1
    imshow(Signal,[])
    title('enhanced signal-image via "top hat" transformation (evens out the background)')
    pause
end
% "Top-hat": substract an opened image from the original. Opening with a
% large enough structural element (disk of radius at least 5 pixels) ensures that
% the opening of the image produces a reasonable estimate of the background


[histy,histx]=imhist(Signal(Signal>0));

%hist(Signal(Signal>0))

[widthx, maxminvalue] = fwhm(histx,histy);
if widthx>0
threshold_cell_auto=maxminvalue+0.8*widthx; %this seems to kinda work for Adam Mig1-GFP data for cell outlines
else
    threshold_cell_auto=maxminvalue;
end
h = fspecial('gaussian');
GaussSignal=imfilter(Signal,h);
%GaussSignal2 = imtophat(GaussSignal,se);
GaussSignal2=GaussSignal;
sigBW = im2bw(GaussSignal2,threshold_cell_auto);

if output==1
   % bar(histx, histy)
   % title(strcat('threshold=',num2str(threshold_cell_auto)))
   imshow(GaussSignal,[])
   title('Gaussian filtered tophat')
    pause
%       imshow(GaussSignal2,[])
%    title('Top-hat Gaussian filtered tophat')
%     pause
    imshow(sigBW,[])
  %  title('thresholded signal-only image')
  title('thresholded top-hat transformed image')
    pause
end



B1 = imopen(sigBW,strel('disk',1)); % morph. opening with disk of radius 1 pixel.
B2 = bwmorph(B1,'fill',1); % fill holes of size one pixel, apply operation one time (1).
%B3 = bwmorph(B2,'open',1); % open image, apply operation one time.
B3=B2;
if output==1
%     imshow(sigBW,[])
%     title('fill holes of size one pixel - new step')
%     pause
    imshow(B1,[])
    title('morph. opening with disk of radius 1 pixel.')
    pause
    imshow(B2,[])
    title('fill holes of size one pixel, apply operation one time (1).')
    pause
%     imshow(B3,[])
%     title('open image, apply operation one time.')
%     pause
end


    B4 = imerode(B3, [1 1; 1 1]); % erode with a 2x2 square shape.
    B5 = imerode(B4, [1 1; 1 1]); % erode with a 2x2 square shape.
    B6 = bwmorph(B5,'shrink',1); % shrinks objects to points, apply operation one time.
    % B6 is a binary image (matrix) with ones at positions of candidate spots and zeros elsewhere.
    result1 = B6;
    if output==1
        imshow(B4,[])
        title('erode with a 2x2 square shape.')
        pause
        imshow(B5,[])
        title('erode with a 2x2 square shape. again')
        pause
        imshow(B6,[])
        title('B6 is a binary image (matrix) with ones at positions of candidate spots and zeros elsewhere.')
        pause
    end
    %% findspots2 style method, better for some data
    
    Bgnd0 = imopen(frame0,se); % Background estimate. Morphological opening (=erosion followed by dilation) with large enough structural element ("se" does not fit entirely within the foreground).
if output==1
    imshow(Bgnd0,[])
    title('Background estimate. Morphological opening (=erosion followed by dilation) with large enough structural element ("se" does not fit entirely within the foreground).')
    pause
end


SignalMask = im2bw(Bgnd0,graythresh(Bgnd0)); % graythresh gives the threshold and im2bw converts image into a thresholded black and white (bw) image.

if output==1
    imshow(SignalMask,[])
    title('thresholded image')
    pause
end
Signal3 = Signal.*SignalMask; % enhanced signal-only image.
if output==1
    imshow(Signal3,[])
    title('enhanced signal-only image')
    pause
end
% Thresholding the signal-only image to try and find spot candidates:
sigBW = im2bw(Signal3,graythresh(Signal3)); % thresholding the top half
if output==1
    imshow(sigBW,[])
    title('thresholded signal-only image')
    pause
end


B1 = imopen(sigBW,strel('disk',1)); % morph. opening with disk of radius 1 pixel.
B2 = bwmorph(B1,'fill',1); % fill holes of size one pixel, apply operation one time (1).
B3 = bwmorph(B2,'open',1); % open image, apply operation one time.
if output==1
    imshow(B1,[])
    title('morph. opening with disk of radius 1 pixel.')
    pause
    imshow(B2,[])
    title('fill holes of size one pixel, apply operation one time (1).')
    pause
    imshow(B3,[])
    title('open image, apply operation one time.')
    pause
end

    B4 = imerode(B3, [1 1; 1 1]); % erode with a 2x2 square shape.
    B5 = imerode(B4, [1 1; 1 1]); % erode with a 2x2 square shape.
    B6 = bwmorph(B5,'shrink',1); % shrinks objects to points, apply operation one time.
    % B6 is a binary image (matrix) with ones at positions of candidate spots and zeros elsewhere.
    result2 = B6;
    if output==1
        imshow(B4,[])
        title('erode with a 2x2 square shape.')
        pause
        imshow(B5,[])
        title('erode with a 2x2 square shape. again')
        pause
        imshow(B6,[])
        title('B6 is a binary image (matrix) with ones at positions of candidate spots and zeros elsewhere.')
        pause
    end

resultTemp=result1+result2;
 resultTemp(resultTemp>1)=1;
 result=bwmorph(resultTemp,'shrink',1);


end