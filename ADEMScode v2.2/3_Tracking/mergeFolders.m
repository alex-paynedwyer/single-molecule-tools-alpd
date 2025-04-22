function mergeFolders(targetfolder,filegrouping)
%run on Sample# subfolder from Micromanager on M1 or M2 to collate acquisitions by field of view
if nargin<1
    targetfolder = uigetdir()
    filegrouping = 2; %brightfield + fluorescence
end
cd(targetfolder)

%filegrouping = 2; %brightfield + fluorescence
%filegrouping = 3; %brightfield + prebleach ch1 + ALEX

makefieldfolder=1;

folders=dir('*_*');
folders = folders(~contains({folders.name},'DS_Store'));
folders = folders(~contains({folders.name},'.ini'));
[folderssort,index]=natsort({folders.name});

% get parent directory (so can be deleted later)
for i = 1:(length(folders))
    % move file to parent directory
    foldername=folderssort{i};
    cd(targetfolder);
    cd(foldername);
    filefolder = dir('*.tif');
    filename = filefolder(1).name;
    metafolder = dir('*metadata.txt');
    commentfolder = dir('*comments.txt');
    dispfolder = dir('*DisplaySettings.json');
    if isempty(metafolder) ~= 1
        metaname = metafolder(1).name;
    end
    if isempty(commentfolder) ~= 1
        commentname = commentfolder(1).name;
    end
    if isempty(dispfolder) ~= 1
        dispname = dispfolder(1).name;
    end

    time=filefolder(1).date;            % rename file to add in time stamp
    time=strrep(time,' ','_');
    time=strrep(time,'-','_');
    time=strrep(time,':','_');
    time=strcat(time(end-7:end),'_',time(1:end-9));
    
    newname=strcat('t_',time,'_',filename);
    if contains(newname,'_MMStack_Pos0')
        newname=replace(newname,'_MMStack_Pos0','');
    elseif contains(newname,'_MMStack.ome')
        newname=replace(newname,'_MMStack.ome','');
    end
    if isempty(metafolder) ~= 1
        newmetaname=strcat('t_',time,'_',metaname);
        if contains(newmetaname,'_MMStack_Pos0')
            newmetaname=replace(newmetaname,'_MMStack_Pos0','');
        elseif contains(newmetaname,'_MMStack')
            newmetaname=replace(newmetaname,'_MMStack','');
        end
    end


if rem(i+filegrouping-1,filegrouping)==0
    cellno = 1+(i-1)/filegrouping;
    if makefieldfolder==1
        newfolder = strcat('../field',num2str(cellno));
        mkdir(newfolder);
    else
        newfolder = '..';
    end
    
end

movefile(filename,strcat(newfolder,'/',newname))   % notify user that empty folder is deleted
if isempty(metafolder) ~= 1
    movefile(metaname,strcat(newfolder,'/',newmetaname))
end
if isempty(commentfolder) ~= 1
    delete(commentname)
end
if isempty(dispfolder) ~= 1
    delete(dispname)
end
disp(['moved file ',filename])
cd ..
rmdir(foldername)  % delete (now) empty parent folder
cd ..
end
end