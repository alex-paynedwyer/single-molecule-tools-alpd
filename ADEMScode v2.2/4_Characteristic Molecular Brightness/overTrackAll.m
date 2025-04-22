%%Overtrack all spots in a sample
% Does not correct for Ch2->Ch1 registration. Only for use in single colour comparisons.

% Directory where data is, can loop through whole strain
dirselect = 1;

if dirselect
    DataDir=uigetdir('*.mat','Select the strain folder to analyse')
    noDataDirs = 1;
else
    DataDirs=["D:\Alex\MATLAB-current\Data\E.coli\...","D:\Alex\MATLAB-current\Data\E.coli\..."];
    noDataDirs = length(DataDirs);
end

AllBaselineCh1=[];
AllBaselineCh2=[];

for i=1:noDataDirs
    if dirselect == 0
        DataDir=DataDirs(i)
    end
    cd(DataDir)
    % Directory to save the analysis data, duplicates hierarchy of data
    AnalysisPath=strcat(DataDir,'_ANALYSIS');
    AltResultsPath=strcat(pwd,'\TrackingFiles');
    DateDir=dir('*2*');          
    DateDir=DateDir([DateDir.isdir]>0);
    %[~, idx] = unique({DateDir.name}.', 'rows', 'stable');  %filter out duplicates
    %DateDir=DateDir(idx);
    %DateDir=DateDir(~ismember({DateDir.name},{'.','..'}));
    DateNum=size(DateDir);

    % Loop over different day folders
    for j=1:DateNum(1) %1:DateNum(1) %change if needed
        cd(strcat(DateDir(j).folder,'\',DateDir(j).name))
        SampleDir=dir('*am*');
        SampleDir=SampleDir([SampleDir.isdir]>0);
        SampleNum=size(SampleDir);
        DatePath=strcat(DataDir, '\', DateDir(j).name);
        DateAPath=strcat(AnalysisPath, '\', DateDir(j).name);
        % Loop over sample folders
        for k=1:SampleNum(1)
            cd(SampleDir(k).name)
            CellDir=dir('*el*');
            CellDir = CellDir(~contains({CellDir.name},'DS_Store'));
            CellDir=CellDir([CellDir.isdir]>0);
            [CellDirSorted,cellindexlist]=natsort({CellDir.name});
            CellNum=size(CellDir);
            SamplePath=strcat(DatePath, '\', SampleDir(k).name);
            SampleAPath=strcat(DateAPath, '\', SampleDir(k).name);
            % Loop over cell folders containing tifs
            for l=1:CellNum(1)    %natural order not alphanumeric
                BFindex=[];
                cd(CellDirSorted{l})
                CellPath=strcat(SamplePath, '\', CellDirSorted{l});
                CellAPath=strcat(SampleAPath, '\', CellDirSorted{l});
                TifFiles=dir('*.tif');
                % Check there are some tifs to analyse
                EnoughTifs=size(TifFiles);
                if EnoughTifs(1)>0
                    %  Check there is a full acquisition by size of file
                    if max([TifFiles(:).bytes]) > 3005
                        pwd
                        FLUORname=TifFiles(find([TifFiles.bytes]== max([TifFiles(:).bytes]))).name;
                        TifFile=strcat(CellPath,'\',FLUORname);
                        DataFile=strcat(CellAPath,'\',FLUORname(1:end-4),'_TRACKS.mat');
                        [spotsBaselineCh1,spotsBaselineCh2]=overTracker(FLUORname,DataFile,0,0,CellPath,CellAPath);
                        AllBaselineCh1=cat(1,AllBaselineCh1,spotsBaselineCh1);
                        AllBaselineCh2=cat(1,AllBaselineCh2,spotsBaselineCh2);
                    end
                end
                cd(CellPath)
                cd ..
            end
            cd ..
        end
        cd ..
    end
end


