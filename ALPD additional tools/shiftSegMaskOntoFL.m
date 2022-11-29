%% CALCULATE & APPLY SHIFT between SEGMENTED MASK and FLUORESCENT FOV

function [shiftedmaskset, bf2fl_tform]=shiftSegMaskOntoFL(maskset,FL)
%maskset ~ CellObject
maskslice = sum(maskset,3); %collapse all slices onto one mask

if size(maskslice)~=size(FL)
    error("Inputs must be the same size.")
end
%Cross correlate the two to find the mean transform vector
c = normxcorr2(FL,maskslice);
[max_c,imax] = max(abs(c(:)));
[ypeak,xpeak] = ind2sub(size(c),imax(1));
tvector = [(xpeak-size(FL,2)) (ypeak-size(FL,1))];

bf2fl_tform = tvector;
shiftedmaskset = circshift(maskset,-wrev(bf2fl_tform));

%bf2fl_tform=affine2d([tvector(1),0,0;0,tvector(2),0;0,0,1]);
%shiftedmaskset = imwarp(maskset,bf2fl_tform,'OutputView',imref2d(size(maskset)));
end