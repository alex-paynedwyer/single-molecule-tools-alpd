
% Loops through a directory containing segmentation data, plots it over brightfield and fluorescence images 
% also fits 'sausage' shapes (ideal for rod-shaped bacteria) and plots those

DefaultPixelSize = 53; %nm
applychanneltransform=0;

% Directory where data is, can loop through whole strain
DataDir=uigetdir('*.mat','Select the strain folder to analyse');
cd(DataDir)

% Directory to save the analysis data, duplicates hierarchy of data
% directory
ResultsPath=strcat(DataDir,'_ANALYSIS\');
%ResultsPath='C:\data\yeast-mig1-nrd1_ANALYSIS\';
AltResultsPath=strcat(pwd,'\TrackingFiles');
%mkdir(AltResultsPath)
%mkdir(ResultsPath)
LoopCount=0;
%directories
DateDir=dir('*2*');
DateNum=size(DateDir);
if DateNum==0
    error('No date folders found.')
end

OffsetPixelNum=13;
pause on
% Loop over different day folders
for j=1:DateNum(1)
    cd(DateDir(j).name)
    SampleDir=dir('*am*');
    SampleNum=size(SampleDir);
    DatePath=strcat(ResultsPath, '\', DateDir(j).name);
    %mkdir(DatePath)
    
    % Loop over sample folders
    for k=1:SampleNum(1)
        cd(SampleDir(k).name)
        CellDir=dir('*el*'); CellDir = CellDir(~contains({CellDir.name},'DS_Store'));
        CellNum=size(CellDir);
        [CellDirSorted,cellindexlist]=natsort({CellDir.name});
        SamplePath=strcat(DatePath, '\', SampleDir(k).name);
        areatot = [];
        numbersegs = 0;
        % mkdir(SamplePath)
        % Loop over cell folders containing tifs
        
        %list4plot=[2];
        list4plot=1:CellNum(1);  %cell folder indices for plot (natural order '1,2,3..', not default alphanumeric order '1,10,11...')
        %list4plot=[1,2,3,4]%,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31];
        z=0;
        for l=1:CellNum(1)  %beware - doesn't include cell0 if exists
            cellindexstr = CellDirSorted{l};
            if contains(cellindexstr,'cell')
                cellindexstr = cellindexstr(5:end);
            elseif contains(cellindexstr,'field')
                cellindexstr = cellindexstr(6:end);
            end
            BFindex=[];
            cellindex = str2num(cellindexstr);
            cd(CellDirSorted{l})
            if sum(ismember(list4plot,cellindex))>0
                z=z+1;
                disp(strcat("Plotting segments for ",CellDirSorted{l}))
                CellPath=strcat(SamplePath, '\', CellDirSorted{l});
                % mkdir(CellPath)
                
                clear xoffset2 yoffset2
                
                TifFiles=dir('*.tif');
                %BF=imread('BF.tif');
                %BF=imread(TifFiles(1).name);
                BF=imread(TifFiles(find([TifFiles.bytes]== min([TifFiles(:).bytes]))).name);
                BF=mat2gray(BF);
                FL=imread(TifFiles(find([TifFiles.bytes]== max([TifFiles(:).bytes]))).name);
                FL=mat2gray(FL);
                
                %calculate the first fluorescent frame
                FLproj=FL(:,:,1);
                
                OrigPath = pwd;
                SegPath = strcat(ResultsPath,DateDir(j).name,'\',SampleDir(k).name,'\',CellDirSorted{l})
                cd(SegPath);
                segfile=dir('*_segmentation.mat');
                load(segfile(end).name);
                cd(OrigPath);
                
                area = sum(sum(CellObject));
                for i = 1:length(area)
                    areatot = [areatot,area(i)];
                end
                areatot = nonzeros(areatot)';  %disregard any fields with no segmented objects
                numbersegs = length(areatot);
                
                figure('units','normalized','outerposition',[0 0 1 1], 'NumberTitle', 'off', 'Name', CellDirSorted{l});  % maximise figure automatically   % show cell identifier in title
                
                % Cell outline on fluorescent image
                
                subplot(2,2,1)
                imshow(mat2gray(Cellframe0))
                hold on
                
                subplot(2,2,2)
                imshow(mat2gray(BF))
                hold on
                
                subplot(2,2,3)
                imshow(mat2gray(FLproj))
                hold on
                
                subplot(2,2,4)
                imshow(mat2gray(FLproj))
                hold on
                
                if area>0
                    for q=1:size(CellObject,3)
                        CellObject2=CellObject(:,:,q);
                        
                        subplot(2,2,1)
                        [row,col]=find(bwperim(CellObject2)==1);
                        if applychanneltransform==1 && exist('tform')
                            [transcol,transrow]=transformPointsForward(tform, col, row);
                            scatter(transcol,transrow,5);
                            disp('Transformed segments used.')
                        else
                            scatter(col,row,5)
                            
                            subplot(2,2,2)
                            scatter(col,row,5)
                            stats=regionprops(CellObject2,'Area','Centroid','MinorAxisLength','MajorAxisLength','Orientation','MinFeretProperties','MaxFeretProperties');
                            Centroid=[stats.Centroid(2),stats.Centroid(1)];
                            
                            FinalImage=SausageFunction(CellObject2,Centroid,stats.MinorAxisLength,stats.MajorAxisLength,90+stats.Orientation);
                            if  stats.MinorAxisLength>5 && stats.MinorAxisLength<25
                                if  stats.MajorAxisLength>30
                                    overlapInt= sum(sum(abs(FinalImage-CellObject2)));
                                    
                                    subplot(2,2,3)
                                    [row2,col2]=find(bwperim(FinalImage)==1);
                                    scatter(col2,row2,5)
                                    
                                    subplot(2,2,4)
                                    scatter(col2,row2,5)
                                    
                                end
                            end
                        end
                    end
                    %pause
                    
                    % clear xoffset sigmax yoffset sigmay, xoffset2 sigmax2 yoffset2 sigmay2, CellObject Cellframe0 ...
                    % CellThresh CellObject2 NucleusObject2 frame0 CellThresh2 NucThresh2 NucObject Nucframe0 NucThresh
                
                end
                disp("Number of segments per field = "+num2str(numbersegs)+" in "+ num2str(z)+"/"+num2str(length(list4plot)));
                disp("Average segment area = "+num2str(mean(areatot*DefaultPixelSize^2/10^6))+" +/- "+ num2str(std(areatot*DefaultPixelSize^2/10^6)));
            end
            cd ..
        end
        cd ..
    end
end
