%%  Function to import volumetric CZI image data and metadata.
% v1.1 - added start and end frame selection
% v1.0 - basic image frame and voxel size import


function [image_data,metadata,frame_Xsize,frame_Ysize,frame_Tsize,frame_Zsize,voxelUnit,voxelSizeX,voxelSizeY,voxelSizeZ] = openCZI_v1_1(fileName,print_metadata,startFrame,endFrame)
if exist('print_metadata')==0
  print_metadata = 0;  
end

data=bfopen(fileName);
grayScale=data{1,1};
metadata = data{1,2};
color_maps = data{1,3};
omeMeta = data{1,4};

if print_metadata == 1
    metadataKeys = metadata.keySet().iterator();
    for i=1:metadata.size()
      key = metadataKeys.nextElement();
      value = metadata.get(key);
      fprintf('%s = %s\n', key, value)
    end
end

frame_Xsize = omeMeta.getPixelsSizeX(0).getValue(); % image width, pixels
frame_Ysize = omeMeta.getPixelsSizeY(0).getValue(); % image height, pixels
frame_Zsize = omeMeta.getPixelsSizeZ(0).getValue(); % number of Z slices
frame_Tsize = length(grayScale)/frame_Zsize;  % number of Z cycles (approx. timepoints)

voxelUnit = char(omeMeta.getPixelsPhysicalSizeX(0).unit().getSymbol()); % returns the default unit type
%voxelSizeX = omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER).doubleValue(); % in µm
%voxelSizeY = omeMeta.getPixelsPhysicalSizeY(0).value(ome.units.UNITS.MICROMETER).doubleValue(); % in µm
%voxelSizeZ = omeMeta.getPixelsPhysicalSizeZ(0).value(ome.units.UNITS.MICROMETER).doubleValue(); % in µm
voxelSizeX = omeMeta.getPixelsPhysicalSizeX(0).value;
voxelSizeY = omeMeta.getPixelsPhysicalSizeY(0).value;
voxelSizeZ = omeMeta.getPixelsPhysicalSizeZ(0).value;
if exist('startFrame') && exist('endFrame') && (startFrame*endFrame > 0)
    image_data=zeros(frame_Ysize,frame_Xsize,endFrame-startFrame+1,frame_Zsize,'uint16');
    for i=1:frame_Zsize
        for j=startFrame:endFrame
            image_data(:,:,j,i)=grayScale{frame_Zsize*(j-1)+i,1};
        end
    end
else
    image_data=zeros(frame_Ysize,frame_Xsize,frame_Tsize,frame_Zsize,'uint16');
    for i=1:frame_Zsize
        for j=1:frame_Tsize
            image_data(:,:,j,i)=grayScale{frame_Zsize*(j-1)+i,1};
        end
    end
end

end