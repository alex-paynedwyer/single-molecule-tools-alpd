%% Plot tracked spots as superresolved overlay on a frame, 
% and export coordinates for ThunderSTORM rendering

% first load segmentation (optional) and track .mat files

localisationprecision = 40; %nm
pixelSize = 53; %nm
singleemitteronly = 1;
singleemitterthreshold = 5;
Isingle = 100;

SNRthresh = 0.5;
firstframe = 1;
endframe = 40;
markerEdgeAlpha = 1;%0.05

scalef = 10;
showframe = 0;

try
    numCells=size(CellObject,3);
catch
    numCells=1;
    CellObject=ones(size(frame_average));
end

scaled_fr_avg = imresize(frame_average,scalef,'nearest');

if showframe
    imshowpair(scaled_fr_avg,zeros(size(frame_average).*scalef)); hold on
else
    imshowpair(zeros(size(frame_average).*scalef),zeros(size(frame_average).*scalef)); hold on
end

allSpots=[];
for i=1:numCells
spots=SpotsCh1(SpotsCh1(SpotsCh1(:,10)>0,11)>SNRthresh,:);
spots=spots(spots(:,9)>firstframe-1,:);
spots=spots(spots(:,9)<endframe+1,:);
if singleemitteronly
    spots=spots(spots(:,5)<(singleemitterthreshold.*Isingle),:);
end
[row,col]=find(bwperim(CellObject(:,:,i))==1);
polyin = convhull(polyshape({row},{col}));
TFin = isinterior(polyin,spots(:,2),spots(:,1));
spots=spots(TFin==1,:);
spots(:,13)=i;
allSpots=[allSpots;spots];
end
x=allSpots(:,1);y=allSpots(:,2);
hold on
scatter(x*scalef,y*scalef,2*scalef*localisationprecision/pixelSize,'filled','Marker','o','MarkerFaceAlpha',markerEdgeAlpha,'MarkerFaceColor','magenta','MarkerEdgeAlpha',markerEdgeAlpha,'MarkerEdgeColor','magenta');

set(gca,'units','pixels') % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels') % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels


T = table(x,y);
%nameout=strcat('allSpots_',num2str(i),'.csv');
nameout=strcat('allSpots.csv');
writetable(T,nameout);




