function alignedSpotsCh2 = applySpotsTform(SpotsCh1,SpotsCh2,tform,framewidth,convert,showOutput)
% Function to apply a predetermined affine transformation to a set of Spot localisations 
% Typically used to match Channel 2 to Channel 1 in dual colour images
alignedSpotsCh2 = SpotsCh2;
if convert ==1 
spotstform=convertspotstform(tform);
[alignedSpotsCh2(:,1), alignedSpotsCh2(:,2)] = transformPointsForward(spotstform, SpotsCh2(:,1), SpotsCh2(:,2));
else
[alignedSpotsCh2(:,1), alignedSpotsCh2(:,2)] = transformPointsForward(tform, SpotsCh2(:,1), SpotsCh2(:,2));
end

if showOutput==1
    figure;
    subplot(1,2,1)
    scatter(SpotsCh1(:,1),SpotsCh1(:,2),'g')
    hold on
    scatter(SpotsCh2(:,1)-framewidth,SpotsCh2(:,2),'r')
    title('misaligned spots')
    
    subplot(1,2,2)
    scatter(SpotsCh1(:,1),SpotsCh1(:,2),'g')
    hold on
    scatter(alignedSpotsCh2(:,1)-framewidth,alignedSpotsCh2(:,2),'r')
    title('aligned spots')
end

