% use with TranslateXY to translate co-ords from channels with different
% X,Y, theta and magnification.

function [xOff, yOff, mag, angle,tform]=ManualAlignBeadsTransform(image1, image2,showoutput)
if nargin<3
    showoutput =1;
end
image1=imadjust(image1,stretchlim(image1,[0 1]),[]);
image2=imadjust(image2,stretchlim(image2,[0 1]),[]);
[moving_out,fixed_out] = cpselect(image2,image1,'Wait', true);
tform = fitgeotrans(moving_out,fixed_out,'NonreflectiveSimilarity');
Jregistered = imwarp(image2,tform,'OutputView',imref2d(size(image1)));
falsecolorOverlay = imfuse(image1,Jregistered);
if showoutput ==1
imshow(falsecolorOverlay,[])
end
u = [0 1];
v = [0 0];
[x, y] = transformPointsForward(tform, u, v);
dx = x(2) - x(1);
dy = y(2) - y(1);
angle = (180/pi) * atan2(dy, dx);
mag = 1 / sqrt(dx^2 + dy^2);
xOff=x(1);
yOff=y(1);
end