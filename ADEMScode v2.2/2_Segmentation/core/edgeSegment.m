function Mask=edgeSegment(frame,method,sizeThresh)

% same as edge segment but creats binary masks rather than label matrices
if nargin<2
    method='Sobel';
end

edgeFrame=edge(frame,method);

se90 = strel('line', 3, 90);
se0 = strel('line', 3, 0);
BWsdil = imdilate(edgeFrame, [se90 se0]);
BWdfill = imfill(BWsdil, 'holes');
BWnobord = imclearborder(BWdfill, 4);
seD = strel('diamond',3);
BWfinal = imerode(BWnobord,seD);
BWfinal2 = imdilate(BWfinal,seD);
bw4 = bwareaopen(BWfinal2, sizeThresh);
%BWfinal = imerode(BWfinal,seD);
cc = bwconncomp(bw4,4);
%segmentation=labelmatrix(cc);
L_cells = labelmatrix(cc);
% figure, imshow(BWfinal), title('segmented image');
    Mask=zeros(size(frame,1),size(frame,2),max(L_cells(:)));
   % clear Mask
    for lbl = 1:max(L_cells(:)); %cycle through all objects in labe matrix... asumume only one nucleus per cell... could fall down here if problem with extended maxima of course!   

            Mask(:,:,lbl) = L_cells == lbl; %# find pixels belonging to current label. First object is the borders I think, so take the objects after these as indivudal cell masks  
    end
    Mask=Mask;
end