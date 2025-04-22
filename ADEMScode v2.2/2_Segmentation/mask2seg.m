%function to convert binary mask from TIFF image file (e.g. ImageJ) to segmentation.mat file for use in track analysis;

function mask2seg()
if nargin<1
    [maskImageFile,maskImageDir]=uigetfile('*.tif');
end
cd(maskImageDir)
[masknumSegs, image_Y, image_X, CellObject, ~] = ExtractImageSequence3(maskImageFile(1:end-4), 1, 1, 1);

CellObject(CellObject>0)=1;

%save empty placeholders for nonexistent brightfield, frame-average etc.
Cellframe0 = zeros(image_Y,image_X);
BF = Cellframe0;
p = [];
p.origin = "Converted directly from TIFF";
bf2fl_tform = affine2d([1,0,0;0,1,0;0,0,1]);  %identity transform

area=squeeze([sum(sum(CellObject,1),2)])';

datafilename=strcat(maskImageFile(1:end-4),'_segmentation.mat');
save(datafilename,'CellObject', 'Cellframe0','BF','p','area','bf2fl_tform')
disp('Segmentation saved from TIFF.')