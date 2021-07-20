function mergeFolders(targetfolder)
%run on Sample# subfolder from Micromanager on M1 or M2 to collate acquisitions by field of view
if nargin<1
    targetfolder = uigetdir()
end
cd(targetfolder)
    
filegrouping = 2; %brightfield + fluorescence
%filegrouping = 3; %brightfield + prebleach ch1 + ALEX

folders=dir;
[folderssort,index]=natsort({folders.name});

% get parent directory (so can be deleted later)
for i = 1:(length(folders)-2)
    % move file to parent directory
    foldername=folderssort{i+2};
    cd(foldername);
    filefolder = dir('*.tif');
    filename = filefolder(1).name;
    metafolder = dir('*metadata.txt');
    if isempty(metafolder) ~= 1
        metaname = metafolder(1).name;
    end
    
    time=filefolder(1).date;            % rename file to add in time stamp
    time=strrep(time,' ','_');
    time=strrep(time,'-','_');
    time=strrep(time,':','_');
    time=strcat(time(end-7:end),'_',time(1:end-9));
    
    newname=strcat('t_',time,'_f',filename);
    if isempty(metafolder) ~= 1 
        newmetaname=strcat('t_',time,'_f',metaname);
    end
    
    if rem(i+filegrouping-1,filegrouping)==0
        cellno = 1+(i-1)/filegrouping;
        newfolder = strcat('../field',num2str(cellno));
        mkdir(newfolder);
    end
    
    movefile(filename,strcat(newfolder,'/',newname))   % notify user that empty folder is deleted
    if isempty(metafolder) ~= 1
        movefile(metaname,strcat(newfolder,'/',newmetaname))
    end
    disp(['moved file ',filename])
    cd ..
    rmdir(foldername)  % delete (now) empty parent folder
end