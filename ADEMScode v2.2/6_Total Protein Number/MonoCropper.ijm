//MonoCropper
//Alex Payne-Dwyer, LEAKE LAB, York. 2023.

//Performs basic batch preprocessing on a video (e.g. SlimVar) of a single ROI in a single (mono)colour channel:
// - background correction (e.g. for detector pattern)
// - windows to a predetermined coarse ROI near the centre of the image
// - refines ROI using basic Otsu thresholding on the first frame to a single, bright convex object per field of view
// - measures morphological parameters of the object
// - crops video to the object

//Input:           [input].tif - an (alphanumeric) list of fluorescence videos, or an alternating list of brightfield and fluorescence videos
//                               (brightfield images are currently ignored)
//(optional)       [darkbg].tif - a dark Slimfield video of background pattern noise, or an average projection thereof.
//                               Must be the same XY size as [input], e.g. taken with the same microscope settings.

//Outputs:    crop_[input].tif - an ROI-cropped video with ROI as overlay (e.g. ready for tracking in ADEMScode/PySTACHIO)
//                 Results.csv - a table of morphological statistics, one row for each refined ROI
//(optional)  mask_[input].tif - a binary image of each refined ROI (e.g. for conversion to _segmentation.mat in ADEMScode)

// - GUI for hyperparameters
#@ File[] (label="Select the file(s) to be analysed", style = "extensions:tif/tiff") fileList
#@ File(style="directory") outputFolder
#@ String(choices={"Brightfield/Fluorescence", "Fluorescence only"}, style="radioButtonHorizontal") fileordermode
#@ String(label="Channel", choices={"Ch1", "Ch2"}, style="radioButtonHorizontal") channelno  //only use Ch2 to segment out the right-hand channel when the input is the full sensor width
#@ String(label="Background subtraction",choices={"On", "Off"}, style="radioButtonHorizontal") bgflag
#@ File[] (label="Select the file to use as dark background", style = "extensions:tif/tiff") bgfileList
#@ Integer(label="Crop x", value=300) cropx //set the location and size of the  initial window to search for a bright object
#@ Integer(label="Crop y", value=128) cropy //
#@ Integer(label="Crop width", value=240) cropw
#@ Integer(label="Crop height", value=260) croph
#@ String(label="Segmentation by intensity",choices={"On", "Off"}, style="radioButtonHorizontal") segflag
#@ String(label="Save binary masks",choices={"On", "Off"}, style="radioButtonHorizontal") savemasks
#@ Integer(label="Minimum segment size (pixels)", value=2000) minsize // only consider objects larger than this
#@ Integer(label="Segment dilation (pixels)", value=20,min=0, max=50, style="slider") segdil // used to set the crop margin around an object
#@ Integer(label="Segment contraction (pixels)", value=20,min=0, max=50, style="slider") segdil
#@ String(label="Save images in separate folders",choices={"On", "Off"}, style="radioButtonHorizontal") saveasSeparateFolders


//#@ String(label="Exp. Group", choices={"Mutant", "Control"}, style="list") expGroup
//#@ String(label="Your Text") userText
//#@ String(value="Some useful hints...", visibility="MESSAGE") hints
//#@ String(label="Analyst name", description="Your name") analystName
//label="A real number") realNumber

//Runtime initialisation of 

if (fileordermode == "Brightfield/Fluorescence") {
start=1;
period = 2;
} else {
start=0;
period = 1;
}

if (channelno == "Ch1") {
xshift = cropx;
} else {
sensorsize = 1200;           //sensorsize in pixels is hardcoded for Photometrics Prime95b
xshift = cropx+sensorsize/2; //expect a second channel to be shifted laterally halfway across the sensor;
}

if (bgflag == "On") {
for (i=0;i<lengthOf(bgfileList);i++){
bgfile = bgfileList[i];
open(bgfile);
}

bgname=getTitle();
if (nSlices>1){
rename("BGstack");
run("Z Project...", "projection=[Average Intensity]");
close("BGstack");
}
rename("BG");
run("Select All");
run("Measure");
bgavg = getResult("Mean",nResults-1);
Table.deleteRows(nResults-1, nResults-1); //delete the last row just measured for dark background
} else {
}

File.makeDirectory(outputFolder+"\\Sample1");
if (savemasks == "On"){
File.makeDirectory(outputFolder+"\\masks");
} else {
}

//Loop over the input file list
file = fileList[0];
j=0;
for (i=start;i<lengthOf(fileList);i=i+period){
file = fileList[i];
j=j+1;
open(file);

imagename=getTitle();
rootname=substring(imagename, 0, lengthOf(imagename)-21);
rename("InputStack");
if (bgflag == "On") {
run("Add...", "value="+bgavg+" stack");
imageCalculator("Subtract create stack","InputStack","BG");
close("InputStack");
selectWindow("Result of InputStack");
rename("InputStack");
} else {
}

run("Set... ", "zoom=100");
setSlice(1);
resetMinAndMax();
Stack.setXUnit("px");
run("Properties...", "pixel_width=1 pixel_height=1 voxel_depth=1 global");
run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret's integrated median skewness kurtosis area_fraction stack display redirect=None decimal=5");
wait(10);
//Crop to specified central region to improve Otsu thresholding of bright feature
run("Specify...", "centered width="+cropw+" height="+croph+" x="+xshift+" y="+cropy+" slice=1");
wait(10);
run("Crop");

if (segflag == "On"){
//setAutoThreshold("Triangle dark");
setAutoThreshold("Otsu dark");
run("Create Mask");
run("Fill Holes");
run("Open");
run("Analyze Particles...", "size="+minsize+"-Infinity pixel add");
wait(10);
selectWindow("mask");
run("Close");
wait(10);
roiManager("Select", 0);
run("Convex Hull");
run("Create Mask");
roiManager("Reset");
run("Create Selection");
wait(10);
roiManager("Add");
selectWindow("Mask");
run("Close");
selectWindow("InputStack");
roiManager("Select", 0);
run("Enlarge...", "enlarge="+segdil+" pixel");
run("Crop");
wait(10);
run("Enlarge...", "enlarge=-"+segdil+" pixel");
//run("Fit Spline");
roiManager("Add");
roiManager("Deselect");
roiManager("Select", 1);
wait(50);
rename(rootname);
roiManager("Measure");
} else {
	rename(rootname);
}

if (savemasks == "On"){
roiManager("Select", 1);
run("Create Mask");
wait(10);
saveAs("TIFF",outputFolder+"\\masks\\mask"+j+"_"+rootname);
rename("MaskOut");
wait(50);
close("MaskOut");
} else {
}
roiManager("Reset");

selectWindow(rootname);
if (saveasSeparateFolders == "On"){
File.makeDirectory(outputFolder+"\\Sample1\\field"+j);
saveAs("TIFF",outputFolder+"\\Sample1\\field"+j+"\\monoC_"+rootname);
} else {
saveAs("TIFF",outputFolder+"\\Sample1\\monoC_"+rootname);
}

selectWindow("ROI Manager");
run("Close");
close();
}
close("BG");
saveAs("Results", outputFolder+"\\Results.csv");
close("Results");