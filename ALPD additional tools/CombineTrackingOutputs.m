%%function to combine TrackAnalyser pairs of output.mat independently of params such as Isingle
function CombineTrackingOutputs(pathA,pathB,saveOutput)

%specify the existing track analysis files and load them
if nargin<1
    saveOutput=1;
end
if nargin<2
    [fileA,pathA] = uigetfile;
    [fileB,pathB] = uigetfile;
end

%if ~exist('pathA') || ~exist('pathB') || ~contains('pathA','output') || ~contains('pathB','output')
%    error('Paths specified must be output.mat files.')
%end

cd(pathA)
load(fileA)
inputA = output;
if ~exist('params')
    params={paramsA,paramsB};
end
params1 = params;
try
    NCellsA = inputA.NCells;
catch
    NCellsA = inputA.Ncells;
end
clear('output','params')

cd(pathB)
load(fileB)
inputB = output;
if ~exist('params')
    params={paramsA,paramsB};
end
params2 = params;
try
    NCellsB = inputB.NCells;
catch
    NCellsB = inputB.Ncells;
end
clear('output','params')

try
    if inputA.NunlinkTrack(2)>0 && inputB.NunlinkTrack(2)>0
        Ch2exists=1;
    else
        Ch2exists=0;
    end
catch
    Ch2exists=0;
end

try
    inputA.AllSpotsA;
    FileAoriginal = 0;
catch
    FileAoriginal = 1;
end
try
    inputB.AllSpotsA;
    FileBoriginal = 0;
catch
    FileBoriginal = 1;
end

%write fields for new combined output

output.NFields = inputA.NFields + inputB.NFields;
output.NCells = NCellsA + NCellsB;
output.NlinkTrack = inputA.NlinkTrack + inputB.NlinkTrack;
output.NunlinkTrack = inputA.NunlinkTrack + inputB.NunlinkTrack;
output.NtotalTrack = output.NlinkTrack + output.NunlinkTrack;
output.linkTrackMean = NCellsA/output.NCells*inputA.linkTrackMean + NCellsB/output.NCells*inputB.linkTrackMean;
output.unlinkTrackMean = NCellsA/output.NCells*inputA.unlinkTrackMean + NCellsB/output.NCells*inputB.unlinkTrackMean;
output.totalTrackMean = output.linkTrackMean + output.unlinkTrackMean;

% container for individual outputs, but not combined properly
if FileAoriginal
    spotsA1=inputA.AllSpots1;
    spotsA2=inputA.AllSpots2;
    tracksA1=inputA.trackArrayCh1;
    tracksA2=inputA.trackArrayCh2;
else
    spotsA1=inputA.AllSpotsA;
    spotsA2=inputA.AllSpotsB;
    tracksA1=inputA.trackArrayA;
    tracksA2=inputA.trackArrayB;
end

if FileBoriginal
    spotsB1=inputB.AllSpots1;
    spotsB2=inputB.AllSpots2;
    tracksB1=inputB.trackArrayCh1;
    tracksB2=inputB.trackArrayCh2;
else
    spotsB1=inputB.AllSpotsA;
    spotsB2=inputB.AllSpotsB;
    tracksB1=inputB.trackArrayA;
    tracksB2=inputB.trackArrayB;
end

output.AllSpotsA = [spotsA1; spotsA2];
output.AllSpotsB = [spotsB1; spotsB2];
output.trackArrayA = [tracksA1; tracksA2];
output.trackArrayB = [tracksB1; tracksB2];
output.alllinkTrack = [inputA.alllinkTrack; inputB.alllinkTrack];
output.allunlinkTrack = [inputA.allunlinkTrack; inputB.allunlinkTrack];

output.segArray = [{inputA.segArray},{inputB.segArray}];

% try combining them properly (e.g. need to add a constant to the trajectory numbers in B to follow on from A)
%output.AllSpots1 = concatenate(inputA.AllSpots1, inputB.AllSpots1);
%output.AllSpots2 = concatenate(inputA.AllSpots2, inputB.AllSpots2);
%output.trackArrayCh1 = concatenate(inputA.trackArrayCh1, inputB.trackArrayCh1);
%output.trackArrayCh1 = concatenate(inputA.trackArrayCh2, inputB.trackArrayCh2);

try  %nest all previous combinations
    output.AllSpotsArray = [inputA.AllSpotsA, inputA.AllSpotsB, inputB.AllSpotsA, inputB.AllSpotsB];
catch
end

output.CellNoofTrajList{1,1} = [inputA.CellNoofTrajList{1,1};inputB.CellNoofTrajList{1,1}];
output.UnlinkedStoichsList{1,1} = [inputA.UnlinkedStoichsList{1,1};inputB.UnlinkedStoichsList{1,1}];
output.LinkedStoichsList{1,1} = [inputA.LinkedStoichsList{1,1};inputB.LinkedStoichsList{1,1}];
output.UnlinkedDiffsList{1,1} = [inputA.UnlinkedDiffsList{1,1};inputB.UnlinkedDiffsList{1,1}];
output.LinkedDiffsList{1,1} = [inputA.LinkedDiffsList{1,1};inputB.LinkedDiffsList{1,1}];
output.UnlinkedTrajList{1,1} = [inputA.UnlinkedTrajList{1,1};inputB.UnlinkedTrajList{1,1}+NCellsA(1)];
output.LinkedTrajList{1,1} = [inputA.LinkedTrajList{1,1};inputB.LinkedTrajList{1,1}+NCellsA(1)];

if Ch2exists
output.CellNoofTrajList{1,2} = [inputA.CellNoofTrajList{1,2};inputB.CellNoofTrajList{1,2}];
output.UnlinkedStoichsList{1,2} = [inputA.UnlinkedStoichsList{1,2};inputB.UnlinkedStoichsList{1,2}];
output.LinkedStoichsList{1,2} = [inputA.LinkedStoichsList{1,2};inputB.LinkedStoichsList{1,2}];
output.UnlinkedDiffsList{1,2} = [inputA.UnlinkedDiffsList{1,2};inputB.UnlinkedDiffsList{1,2}];
output.LinkedDiffsList{1,2} = [inputA.LinkedDiffsList{1,2};inputB.LinkedDiffsList{1,2}];
output.UnlinkedTrajList{1,2} = [inputA.UnlinkedTrajList{1,2};inputB.UnlinkedTrajList{1,2}+NCellsA(2)];
output.LinkedTrajList{1,2} = [inputA.LinkedTrajList{1,2};inputB.LinkedTrajList{1,2}+NCellsA(2)];
output.separationMean(1) = (NCellsA/output.NCells*inputA.separationMean(1)+NCellsB/output.NCells*inputB.separationMean(1)); 
output.separationMean(2) = (NCellsA/output.NCells*inputA.separationMean(1)+NCellsB/output.NCells*inputB.separationMean(1))*sqrt((NCellsA/output.NCells)*inputA.separationMean(2)^2/inputA.separationMean(1)^2+(NCellsB/output.NCells)*inputB.separationMean(2)^2/inputB.separationMean(1)^2);

if exist('inputA.PairedStoichsList')&&exist('inputB.PairedStoichsList')
output.PairedStoichsList = [inputA.PairedStoichsList;inputB.PairedStoichsList];
output.PairedDiffsList = [inputA.PairedDiffsList;inputB.PairedDiffsList];
end

output.unlinkStoichsMean = [mean(output.UnlinkedStoichsList{1,1}), mean(output.UnlinkedStoichsList{1,1})/(output.NunlinkTrack(1)^0.5);...
    mean(output.UnlinkedStoichsList{1,2}), mean(output.UnlinkedStoichsList{1,2})/(output.NunlinkTrack(2)^0.5)];
output.linkStoichsMean = [mean(output.LinkedStoichsList{1,1}), mean(output.LinkedStoichsList{1,1})/(output.NlinkTrack(1)^0.5);...
    mean(output.LinkedStoichsList{1,2}), mean(output.LinkedStoichsList{1,2})/(output.NlinkTrack(2)^0.5)];
output.unlinkDiffsMean = [mean(output.UnlinkedDiffsList{1,1}), mean(output.UnlinkedDiffsList{1,1})/(output.NunlinkTrack(1)^0.5);...
    mean(output.UnlinkedDiffsList{1,2}), mean(output.UnlinkedDiffsList{1,2})/(output.NunlinkTrack(2)^0.5)];
output.linkDiffsMean = [mean(output.LinkedDiffsList{1,1}), mean(output.LinkedDiffsList{1,1})/(output.NlinkTrack(1)^0.5);...
    mean(output.LinkedDiffsList{1,2}), mean(output.LinkedDiffsList{1,2})/(output.NlinkTrack(2)^0.5)];

else  %throw away second channel if any Ch2 is empty

output.unlinkStoichsMean = [mean(output.UnlinkedStoichsList{1,1}), mean(output.UnlinkedStoichsList{1,1})/(output.NunlinkTrack(1)^0.5)];
output.linkStoichsMean = [mean(output.LinkedStoichsList{1,1}), mean(output.LinkedStoichsList{1,1})/(output.NlinkTrack(1)^0.5)];
output.unlinkDiffsMean = [mean(output.UnlinkedDiffsList{1,1}), mean(output.UnlinkedDiffsList{1,1})/(output.NunlinkTrack(1)^0.5)];
output.linkDiffsMean = [mean(output.LinkedDiffsList{1,1}), mean(output.LinkedDiffsList{1,1})/(output.NlinkTrack(1)^0.5)];

end 

paramsA = params1;
paramsB = params2;

if saveOutput==1
        cd(pathA)
        savefile = 'combined_output.mat'; %simple name
        save(savefile,'output','paramsA','paramsB');
end

end
