%% PREPARE IMAGE REGISTRATION for spotChSimilarity based on shift in brightfield

function [strregframech1,regframech2,framech1,framech2,CellObject]=prepRegImages(BFImage,p)

if ~exist('p.CSplit')
    p.CSplit=1;
end
if ~exist('p.bitDepth')
    p.bitDepth=12;
end
if ~exist('p.topbottomcut') || ~exist('p.leftrightcut')
topbottomcut = 5;%5;
leftrightcut = 70;%70;
end
gaussweight=1;
    
unalignedImage = BFImage;
if p.CSplit==2
    unalignedImage = rot90(unalignedImage,3);  %rotation necessary for horizontal camera slits
end
[image_Y,image_X,numFrames]=size(unalignedImage);
frameheight=round(image_Y);
framewidth=round(image_X/2);
framech1=unalignedImage(:,1:framewidth,1);
framech2=unalignedImage(:,1+framewidth:end,1);
%figure; subplot(1,3,1); imshow(framech1,[]); subplot(1,3,2); imshow(framech2,[]);subplot(1,3,3);imshowpair(framech1,framech2);

%eliminate black slit areas outside channels from registration weighting by first blurring the whole image
if max(topbottomcut,leftrightcut)>0
cropframech1=imgaussfilt(framech1,max(topbottomcut,leftrightcut));
cropframech2=imgaussfilt(framech2,max(topbottomcut,leftrightcut));

%replace the sharply resolved central areas so these determine the registration
cropframech1(1+topbottomcut:end-topbottomcut,1+leftrightcut:end-leftrightcut)=framech1(1+topbottomcut:end-topbottomcut,1+leftrightcut:end-leftrightcut);
cropframech2(1+topbottomcut:end-topbottomcut,1+leftrightcut:end-leftrightcut)=framech2(1+topbottomcut:end-topbottomcut,1+leftrightcut:end-leftrightcut);
else
    cropframech1=framech1;
    cropframech2=framech2;
end

if gaussweight==1
gausswidth = framewidth/2;
gaussstrength = 0.5; %1 = edges are black, 0 = no correction;
[Xreg,Yreg] = meshgrid(1:framewidth,1:frameheight);    
gaussfactor = double(exp(-(((double(Xreg)-framewidth/2).^2)+(double(Yreg)-frameheight/2).^2)./(2.*gausswidth^2)));
strregframech1 = uint16(cropframech1.*(gaussstrength*gaussfactor+1.0-gaussstrength));
regframech2 = uint16(cropframech2.*(gaussstrength*gaussfactor+1.0-gaussstrength));
else
   strregframech1=cropframech1;
   strregframech1=cropframech2;
end

%GENERATE SEGMENTATION PLACEHOLDER
%areach1=SimpleSegment(framech1,1,0.7*mean(framech1,'all')/2^p.bitDepth); %get mask for channel 1 slit region
%areach2=SimpleSegment(framech2,1,0.7*mean(framech2,'all')/2^p.bitDepth); %get mask for channel 2 slit region
%CellObject=[sum(areach1,3),sum(areach2,3)]; %save channel slit regions as guess for unregistered segmentation window
CellObject=ones(size(BFImage));
%figure; imshow(CellObject(:,:,1),[])

framech1=histeq(uint16(framech1));
framech2=histeq(uint16(framech2));

strregframech1=histeq(uint16(strregframech1));
regframech2=histeq(uint16(regframech2));
end
