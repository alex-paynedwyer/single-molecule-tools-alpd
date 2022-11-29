function placeIntoFolders(targetfolder)
%run on Sample# subfolder from Micromanager on M1 or M2 to collate acquisitions by field of view
if nargin<1
    targetfolder = uigetdir()
end
cd(targetfolder)

separate_samples=1; %place into separate SampleX folders instead of fields in Sample1

filegrouping = 1; %fluorescence
%filegrouping = 2; %brightfield + fluorescence
%filegrouping = 3; %brightfield + prebleach ch1 + ALEX

tiffiles=dir('*.tif');
[tiffilessort,index]=natsort({tiffiles.name});

% get parent directory (so can be deleted later)
for i = 1:(length(tiffilessort))
    % move file to parent directory
    filename=tiffilessort{i};
    if rem(i+filegrouping-1,filegrouping)==0
        fovno = 1+(i-1)/filegrouping;
        if separate_samples==1
        newfolder = strcat('./22MMDD/Sample',num2str(fovno),'/field1');
        else
        newfolder = strcat('./22MMDD/Sample1/field',num2str(fovno));
        end
        mkdir(newfolder);
    end
    movefile(filename,strcat(newfolder,'/',filename))
    matfiles=dir(strcat('*',filename(1:end-4),'*_GRNDTRUTH.mat'));
    if isempty(matfiles) ~= 1
        matfilename=matfiles(1).name;
        movefile(matfilename,strcat(newfolder,'/',matfilename))
    end
    disp(['moved file ',filename])
end