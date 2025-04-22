% Function to refine a non-linear transformation between 2 channels of a
% microscope based on similarity registration transforms of e.g. brightfield images,
% starting from an initial transform derived from e.g. manual pairing from bead images.

function [imageCh2start,imageCh2adj,tformout]=transformRefinement(imageCh1, imageCh2, tformin, invertflag, showOutput)

%invertflag not used yet, placeholder if needed for an alternative comparison of ch1 -> ch2 whose output transform could then be inverted; 
% that route may have better optimisation performance depending on the relative contrast in the two channels

imageCh1=histeq(uint16(imageCh1));
imageCh2=histeq(uint16(imageCh2));

[optimizer, metric] = imregconfig('monomodal');

optimizer = registration.optimizer.OnePlusOneEvolutionary;
optimizer.InitialRadius = 0.0001; %0.0001;
optimizer.Epsilon = 1.5e-6; %1.5e-6;
optimizer.GrowthFactor = 1.0005; %1.001;
optimizer.MaximumIterations = 1000; %1000;

metric = registration.metric.MattesMutualInformation;
metric.NumberOfSpatialSamples = 500;    %500;
metric.NumberOfHistogramBins = 100;     %100;
metric.UseAllPixels = 1;

transformtype ='similarity';  %'translation';'rigid';'similarity';'affine';

tformout = imregtform(imageCh2,imageCh1,transformtype,optimizer,metric,'InitialTransformation',tformin); 

% check transform works on dual colour images
imageCh2start = imwarp(imageCh2,tformin,'OutputView',imref2d(size(imageCh1)));
imageCh2adj = imwarp(imageCh2,tformout,'OutputView',imref2d(size(imageCh1)));

if showOutput==1
figure('Name','Original','NumberTitle','off')
falsecolorOverlay1 = imfuse(imageCh1,imageCh2);
imshow(falsecolorOverlay1,[])
figure('Name','Input transformation','NumberTitle','off')
falsecolorOverlay2 = imfuse(imageCh1,imageCh2start);
imshow(falsecolorOverlay2,[])
figure('Name','Output transformation','NumberTitle','off')
falsecolorOverlay3 = imfuse(imageCh1,imageCh2adj);
imshow(falsecolorOverlay3,[])
end

end