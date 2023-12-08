#@File (label="Input folder", style="directory") userChosenDirectory
#@File(label = "Output directory", style = "directory") userChosenOutputDirectory
#@String (choices={"Ignore", "DAPI", "O4", "MBP"}, style="list") channel_1
#@String (choices={"Ignore", "DAPI", "O4", "MBP"}, style="list") channel_2
#@String (choices={"Ignore", "DAPI", "O4", "MBP"}, style="list") referenceChannelTag
#@Integer (label="Min. overlap %", value=5) minOverlapPercentage
#@Boolean (label="Measure morphology", value=false) morphologyMeasureBoolean
#@Boolean (label="Limit image number", value=false) limitSeriesCount
#@Integer (label="Max. image / folder", value=5) maxSeriesNumber
#@String(label = "File suffix", value = ".czi") fileSuffix

//setup
userChosenDirectory += File.separator;
setOption("ScaleConversions", true);
setOption("BlackBackground", true);
run("Bio-Formats Macro Extensions");


//global variables
var seriesCount=0;
var originalImageName;
var shortOriginalImageName;
var extractedImageName;
var outputDirWithOrigImageNameWithStaining;
var activeChannelCounter = 0;
//var referenceChannelTag = "DAPI";
var referenceChannelName = "";
var referenceChannelID = 0;
var referenceChannelIndex;
var referenceObjectCount = 0;
var headings = "Image,Nuclei,";
var dataArray;
var dataTable;
var stainingArray = newArray(0);
var stainingArea;
var	totalBranchLengthSkel;
var totalBranchLengthRD;
var totalJunctionsSkel;
var totalJunctionsRD;
var totalTriplePointsSkel;
var	totalQuadruplePointsSkel;
var directoryHasSubfolder = false;


//Program starts here
processBioFormatFiles(userChosenDirectory);
close("BFE_Results");
close(dataTable);
close("ROI Manager");
exit("Done");
//Execution stops here


//Loop through folder structure
function processBioFormatFiles(currentDirectory) {
	
	continueProcessing = true;
	//currentDirectory = userChosenDirectory += File.separator; so path with a \ at the end. 
	directoryWithoutFileSeparator = substring(currentDirectory, 0, lengthOf(currentDirectory)-1);
	//lastDirectoryName = substring(directoryWithoutFileSeparator, lastIndexOf(directoryWithoutFileSeparator, "_")+1);
	lastDirectoryName = substring(directoryWithoutFileSeparator, lastIndexOf(directoryWithoutFileSeparator, File.separator)+1);
	print("currentDirectory: "+ currentDirectory);
	print("directoryWithFileSeperator: "+substring(currentDirectory, 0, lengthOf(currentDirectory)));
	print("directoryWithoutFileSeparator: " + directoryWithoutFileSeparator);
	print("lastDirectoryName: "+lastDirectoryName);
	print("lasIndexOf: " + lastIndexOf(directoryWithoutFileSeparator, File.separator));
	fileList = getFileList(currentDirectory);
	Array.print(fileList);
	print("fileLIST ");
	//if limitSeriesCount = true, series = images/folder
	if (limitSeriesCount == false) {
			maxSeriesNumber = true; 
	}	
	//loop through all files, or chosen limited amount of images.
	for (file = 0; file < fileList.length && seriesCount < maxSeriesNumber; file++) {
		if (endsWith(fileList[file], fileSuffix)) {
			//increase seriesCount for each image edited, if limit is activated (limitSeriesCount=true)
			if (limitSeriesCount) {
			seriesCount++;
			}
			//if next file doesnt have the same beginning. like B1_1 and B4_1 or B2_3, close dataTable
			//create a new table for the first file and first file of each "series" like B1, B2...; in both cases dataTable is not open
			if(!isOpen(dataTable)){
				//print("New Table Substring file: " + substring(fileList[file], 0, indexOf(fileList[file], "_")));	
				dataTable = "" + getDateStamp() + "_Analysis-" + lastDirectoryName;
				//dataTable = "" + "_Analysis-" + lastDirectoryName;
				Table.create(dataTable);
			}
			Ext.setId(currentDirectory + fileList[file]);
			run("Bio-Formats Importer", "open=[" + currentDirectory + fileList[file] + "] autoscale color_mode=Colorized rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
			processImage();
			saveData();
			close("*");	
			//either for every image, or for limited number of images, do processImage
			//for (series = 1; series <= maxSeriesNumber; series++) {
				//Open Image of series (problem series can be 3 different resolutions of the SAME image
				//run("Bio-Formats Importer", "open=[" + currentDirectory + fileList[file] + "] autoscale color_mode=Colorized rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_" + series);
				//Open Image without series 

			//}
			//if next file doesnt have the same beginning. like B1_1 and B4_1 or B2_3, close dataTable
			if(file>=fileList.length){
				if(seriesCount >= maxSeriesNumber && substring(fileList[file], 0, indexOf(fileList[file], "_")) != substring(fileList[file+1], 0, indexOf(fileList[file+1], "_"))){
					print("Substring file: " + substring(fileList[file], indexOf(fileList[file], "_")));
					print("Substring file+1: "+substring(fileList[file+1], indexOf(fileList[file+1], "_")));
					close(dataTable);
				}
			}
				
		
		} else if (File.isDirectory(currentDirectory + fileList[file])) {
			//if current file is a folder, pass on as new current Directory to be processed
			directoryHasSubfolder = true;
			processBioFormatFiles(currentDirectory + fileList[file]);
		}
	} 
}




function processImage() {
	originalImageID = getImageID();
	//originalImageName e.g.: B1_1.czi - B1_1.czi #1
	originalImageName = getTitle();
	//shortOriginalImageName e.g.: B1_1
	shortOriginalImageName = substring(originalImageName, 0, indexOf(originalImageName, "."));
	//if directory has subfolders, original image name is B1/B1_1
	if(directoryHasSubfolder){
		shortOriginalImageName = substring(originalImageName, indexOf(originalImageName, "/")+1, indexOf(originalImageName, "."));
	}
	run("8-bit");
	getDimensions(width, height, channels, slices, frames);
	dataArray = newArray(8);
	
	//set up empty variables for saving dimensions of channels
	stainingDefinitions = newArray(channel_1, channel_2);
	channelID = newArray(channels);
	channelName = newArray(channels);
	
	activeChannelCounter = 0;
	ignoredChannelCount = 0;
	for (c = 0; c < channels; c++) {
		//check for all channels, whether chosen channel is defined. Then process it. Otherwise ("ignore") delete the channel variable from array
		if (stainingDefinitions[c-ignoredChannelCount] != "Ignore") {
			
			activeChannelCounter++;

			selectImage(originalImageID);
			run("Duplicate...", "duplicate channels=" + (c+1));	
			channelID[c-ignoredChannelCount] = getImageID();
			channelName[c-ignoredChannelCount] = getTitle();
			//if chosen channel == DAPI, set references to properties of this channel
			if (matches(stainingDefinitions[c-ignoredChannelCount], referenceChannelTag)) {
				referenceChannelID = channelID[c-ignoredChannelCount];
				referenceChannelName = channelName[c-ignoredChannelCount];
				referenceChannelIndex = c-ignoredChannelCount;
			}
			//extractObject will edit channel
			extractObject(stainingDefinitions[c-ignoredChannelCount]);
			//if channel doesnt exist, delete it and +1 for variable 
		} else {
			stainingDefinitions = Array.deleteIndex(stainingDefinitions, c-ignoredChannelCount);
			ignoredChannelCount++;
		}
		
	}
	//check whether Reference channel was defined
	if (referenceChannelID == 0 || referenceChannelName =="" || referenceChannelTag =="Ignore" ) {
		exit("Reference staining missing!");
	}

	//Start with processing referenceChannel e.g. DAPI
	selectImage(referenceChannelID);
	selectWindow(referenceChannelName);
	// ANALYZE > SET MEASUREMENTS > redirect to none
	run("Set Measurements...", "area mean redirect=None decimal=3");
	//Count amount of e.g. DAPI-nuclei
	run("Analyze Particles...", "size=30.00-infinity circularity=0.00-1.00 show=Nothing clear add");
	referenceObjectCount = nResults;
	//dataArray[0] = referenceObjectCount;
	//delete pixles outside of ROIs in DAPI image (later binary feature extractor uses this to calculate co-localization
	roiManager("Show All");
	roiManager("Combine");
	run("Clear Outside");
	//if no nuclei detected in roiManager, clear out whole image
	if (referenceObjectCount == 0) {
		run("Select All");
		setBackgroundColor(0, 0, 0);
		run("Clear", "slice");
	}

	Table.reset("BFE_Results");
	
	//single positives, if picture includes only one channel extra to the reference channel
	singlePositiveChannelName = newArray(activeChannelCounter-1);
	singlePositiveChannelTag = newArray(activeChannelCounter-1);
	singlePositiveCounter = 0;
	for (comp1 = 0; comp1 < activeChannelCounter; comp1++) {
		//for all channels but Reference channel count colocalization
		if (comp1 != referenceChannelIndex) {
			//draw circles around DAPI nuclei, in case 2nd staining is MBP (since before subtraction MBP-DAPI)
			//if(stainingDefinitions[comp1]=="MBP"){
			//redefined since MBP strategy has changed
			if(stainingDefinitions[comp1]=="O4"){
				run("Binary Feature Extractor", "objects=["+referenceChannelName+"] selector=["+channelName[comp1]+"] object_overlap="+minOverlapPercentage+" count");
				//rename Extracted image created from binary feautre extractor
				extractedImageName = (originalImageName + "_" + stainingDefinitions[referenceChannelIndex] + "-" + stainingDefinitions[comp1]);
				rename(originalImageName + "_" + stainingDefinitions[referenceChannelIndex] + "-" + stainingDefinitions[comp1]);
				singlePositiveChannelName[singlePositiveCounter] = getTitle();
				saveAs("png", outputDirWithOrigImageNameWithStaining + "_DAPI_positive_nuclei" +".png");
			//If 2nd/AB channel not MBP, but O4, use Binary Feature Extractor
			} else if(stainingDefinitions[comp1]=="MBP") {
				run("Binary Feature Extractor", "objects=["+referenceChannelName+"] selector=[Result of "+stainingDefinitions[comp1]+"_black_white] object_overlap="+minOverlapPercentage+" count");
				//rename Extracted image created from binary feautre extractor
				extractedImageName = (originalImageName + "_" + stainingDefinitions[referenceChannelIndex] + "-" + stainingDefinitions[comp1]);
				rename(originalImageName + "_" + stainingDefinitions[referenceChannelIndex] + "-" + stainingDefinitions[comp1]);
				singlePositiveChannelName[singlePositiveCounter] = getTitle();
				saveAs("png", outputDirWithOrigImageNameWithStaining + "_DAPI_positive_nuclei" +".png");
			} else{
				//Binary Feature Extractor creates a new image including only the particles of the "objects image" that overlap with the "selector image" to a certain percentage (object_overlap)
				//"count" creates a table ("BFE results") with n(Objects), n(Selectors), n(Extracted)
				//analysis creates a table of each image (object, selector and extracted image) with 
				run("Binary Feature Extractor", "objects=["+referenceChannelName+"] selector=["+channelName[comp1]+"] object_overlap="+minOverlapPercentage+" count");
				//print("First value: " + Table.get("Image", 0, "BFE_Results") + "length of table: " + getValue("results.count"));
				
				//first row sometimes contains a false line with zeros. This should be deleted
				for(i = 0; i < Table.size("BFE_Results"); i++){
					if (Table.get("Image", i) == 0){
						Table.get("Image", 1);
						Table.deleteRows(i, i);
						print("first row = 0");
					}
				}
				//rename Extracted image created from binary feautre extractor
				extractedImageName = (originalImageName + "_" + stainingDefinitions[referenceChannelIndex] + "-" + stainingDefinitions[comp1]);
				rename(originalImageName + "_" + stainingDefinitions[referenceChannelIndex] + "-" + stainingDefinitions[comp1]);
				singlePositiveChannelName[singlePositiveCounter] = getTitle();
				saveAs("png", outputDirWithOrigImageNameWithStaining + "_DAPI_positive_nuclei" +".png");
				
			}
			singlePositiveChannelTag[singlePositiveCounter] = stainingDefinitions[comp1];
			//headings = headings in final result table (starts with Image, Nuclei, )
			//add channel name to headings, later used for result table to display amount of colocalized cells.
			//headings += singlePositiveChannelTag[singlePositiveCounter] + ",";
			singlePositiveCounter++;
		}
	}
}


function extractObject(staining) {
	//create variable with path and image name, to save files later
	outputDirWithOrigImageNameWithStaining = userChosenOutputDirectory + "/" + 
	//substring creates image name without extensions like .czi, .fsv
	shortOriginalImageName + "_" + staining;
	var staining_channelID = getImageID();
	var staining_channelName = getTitle();
	
	if (staining == "DAPI") {
		extractDAPI();
	} else if (staining == "O4") {
		extractO4();
	} else if (staining == "MBP") {
		extractMBP();
	}
}


function extractDAPI() {
	run("Pseudo flat field correction", "blurring=1500 hide");
	run("Median...", "radius=2");
	run("Subtract Background...", "rolling=60 sliding");
	run("Enhance Contrast...", "saturated=0.01 normalize");
	run("Auto Threshold", "method=Li white");
	//remove small particles
	run("EDM Binary Operations", "iterations=6 operation=open"); 
	//close holes in nuclei
	run("EDM Binary Operations", "iterations=2 operation=close");	
	run("Watershed");
	
	
	//delete particles that are smaller than a nuclei and dont seem to have a form like a nuclei (solidity <.50)
	run("Set Measurements...", "area redirect=[referenceChannelName] decimal=3");
	run("Extended Particle Analyzer", "area=30-Infinity solidity=0.5-1");
	close("Exception");
	roiManager("Show All");
	roiManager("Combine");
	run("Clear Outside");
	
	//duplicate image to produce ROIs from black and white image
	run("Duplicate...", " ");
	//change image name
	rename("DAPI_channel_for_MBP_signal_extraction");
	
	//increase DAPI signal for later subtraction of MBP signal
	//not needed since new method of mbp editing (
	//run("EDM Binary Operations", "iterations=15 operation=dilate");
	//switch back to original image
	selectWindow(referenceChannelName);
	
	saveExtractedImage();
}

function extractO4(){
	//1. correct uneven lighting (flat-field)
	//The radius should be around a tenth or a twentieth of the size of the image. 
	//For example if the image is 512 by 512 pixels, the blur radius would be 25-50. 
	run("Pseudo flat field correction", "blurring=400 hide");
	//2. Image filtering
	//Filter
	run("Median...", "radius=4");
	//3. Reduce Background: 
	//Process > Substract Background > Sliding paraboloid, radius = diameter of smallest non-background object
	run("Subtract Background...", "rolling=500");
	//4. enhance contrast
	//PROCESS > Enhance Contrast > 
	run("Enhance Contrast...", "saturated=0.1 normalize");
	
	//duplicate image to produce ROIs from black and white image
	run("Duplicate...", " ");
	//change image name
	rename(staining+"_black_white");
	
	//5. Set Tresholding	
	// IMAGE > ADJUST > Treshold
	// Set Treshold to Mean for Antibody(AB) channel (oder Li  für O4)
	//setAutoThreshold("Li dark");	//auto thresholds from the manual thresholding dialog are sometimes different and might be suboptimal
	run("Auto Threshold", "method=Huang ignore_black ignore_white white");
	//6. Apply necessary binary operations
	//EDM binary - Remove small particles
	run("EDM Binary Operations", "iterations=2 operation=open");
	//thin shape
	run("EDM Binary Operations", "iterations=2 operation=erode");
	
	//save copy of image of O4 staining
	saveExtractedImage();
	
	//measure area and save as stainingArea
	run("Set Measurements...", "area redirect=["+staining+"_black_white] decimal=3");
	run("Analyze Particles...", "size=40-Infinity solidity=0.00-1 display clear include summarize add");
	
	//duplicate obtained image for later analysis
	run("Duplicate...", " ");
	rename(staining+"_black_white_2");
	selectWindow(staining+"_black_white");
	
	//measure area and save as O4 Area, exclude artefacts with cirularity<0.8
	run("Extended Particle Analyzer", "  area=40-Infinity solidity=0.00-1 show=Nothing redirect=None keep=None summarize add");
	Table.rename("Summary", staining + "_Area");
	//close created image, (mask of mask of mbp_black_white)
	//close();
	//delete particles that didnt fulfill criteria of MBP signal
	//delete pixles outside of ROIs in black and white image
	roiManager("Show All");
	roiManager("Combine");
	run("Clear Outside");
	
	
	//save size of O4-Area
	stainingArea = Table.getColumn("Total Area", staining+"_Area");
	//close "summary" table
	close(staining+"_Area");
		
	//always add staining name, percent of colocalization and stainingArea to headings if activated, later used to display result table 
	headings += staining + "," + staining + "_per" + "," + staining+"Area_micron^2"+ ",";
	stainingArray = Array.concat(stainingArray, stainingArea);
	
	
	//***************Morphology Analysis**********************
	//calculate morphology of cells only if activated by user
	//delete pixles outside of ROIs in original result image
	selectWindow(staining_channelName);
	run("Subtract Background...", "rolling=50");
	run("Enhance Contrast...", "saturated=0.5");
	
	roiManager("Show All");
	roiManager("Combine");
	run("Clear Outside");
	
	
	rename("Result of " + staining + "_black_white");
	
	if(morphologyMeasureBoolean){
		measureMorphologyFunc();
	}
	//close not edited original MBP image, so that co-lokalization doesnt get confused. (otherwise two images with same name)
	close(staining_channelName);
	//rename O4_black_white to original image name to use it in co-lokalization later
	selectWindow(staining+"_black_white_2");
	
	//selectWindow("Result of " + staining + "_black_white");
	rename(staining_channelName);
	
	//Ridge Detection to capture edge of O4-islands
	//run("Ridge Detection", "line_width=3.5 high_contrast=230 ...
	//low_contrast=87 extend_line displayresults method_for_overlap_resolution=NONE ...
	//sigma=1.51 lower_threshold=1 upper_threshold=2 minimum_line_length=0 maximum=0");
}

function extractMBP(){	
	//duplicate image for later use in ridge detection
	run("Duplicate...", " ");
	rename(staining+"_black_white");
	run("Pseudo flat field correction", "blurring=1500 hide");
	run("Median...", "radius=2");
	
	//run("Subtract Background...", "rolling=50 sliding"); Subtraction is altering the fine MBP signal too much
	//run("Auto Threshold", "method=[Try all] ignore_black ignore_white white");
	run("Auto Threshold", "method=Huang2 ignore_black ignore_white white");
	
	//Measure 
	run("Set Measurements...", "redirect=["+staining_channelName+"] decimal=3");
	//Set up ROIs in MBP image
	run("Extended Particle Analyzer", "area=65-Infinity solidity=0.00-0.7 summarize add");
	close("Exception");
	//delete pixles outside of ROIs in black and white image that dont meet the criteria for being mbp signal
	roiManager("Show All");
	roiManager("Combine");
	run("Clear Outside");
	
	
	//Measure mean gray value in MBP image
	run("Set Measurements...", "area mean redirect=["+staining_channelName+"] decimal=3");
	run("Measure");
	//close Results
	close("Results");
	//check whether valid MBP signal was found
	if(roiManager("count")==0){
		print("roiManager empty");
		run("Select All");
		setBackgroundColor(0, 0, 0);
		run("Clear", "slice");
	} else{
		//Measure each ROI individually
		roiManager("Show All");
		roiManager("Combine");
		roiManager("measure");
		//rename Result table
		Table.rename("Results", "ROI in MBP results");
	}
	

	//Measure mean gray value in ROI areas in DAPI image
	run("Set Measurements...", "area mean redirect=["+referenceChannelName+"] decimal=3");
	run("Measure");
	//close Results
	close("Results");
	//check whether valid DAPI signal was found
	if(roiManager("count")==0){
		print("roiManager empty");
		run("Select All");
		setBackgroundColor(0, 0, 0);
		run("Clear", "slice");
	} else{
		//Measure each ROI individually
		roiManager("Show All");
		roiManager("Combine");
		roiManager("measure");
		//rename Result table to 1st results
		Table.rename("Results", "ROI in DAPI results");	
	}
	
	
	
	
	//Binary Feature Extractor creates a new image including only the particles of the "objects image" that overlap with the "selector image" to a certain percentage (object_overlap)
	//"count" creates a table ("BFE results") with n(Objects), n(Selectors), n(Extracted)
	//analysis creates a table of each image (object, selector and extracted image) with 
	//create image that contains !wrong! MBP Rois defined by overlap of DAPI more than 60%, then subtract extracted MBP image from original MBP_black_white
	//run("Binary Feature Extractor", "objects=["+referenceChannelName+"] selector=["+channelName[comp1]+"] object_overlap="+(100-minOverlapPercentage)+" count");
	run("Binary Feature Extractor", "objects=[MBP_black_white] selector=[DAPI_channel_for_MBP_signal_extraction] object_overlap=60 count");
	imageCalculator("Subtract create", "MBP_black_white", "Extracted_MBP_black_white-1");
	//close extracted image from binary feature extractor (the image that contains false positive MBP signal)
	close("Extracted_MBP_black_white-1");
	close("ROI in MBP results");
	close("ROI in DAPI results");
	
	selectWindow("Result of " + staining + "_black_white");
	//save copy of image
	saveExtractedImage();
	
	//obtain for ROIs
	run("Analyze Particles...", "clear include add");
	
	//Measure mean gray value in ROI areas in DAPI image
	run("Set Measurements...", "area redirect=[Result of " + staining + "_black_white] decimal=3");
	run("Measure");
	Table.rename("Summary", staining + "_Area");

	//save size of MBP-Area
	stainingArea = Table.getColumn("Total Area", staining+"_Area");
	//close "summary" table
	close(staining+"_Area");
		
	//always add staining name, percent of colocalization and stainingArea to headings if activated, later used to display result table 
	headings += staining + "," + staining + "_per" + "," + staining+"Area_micron^2"+ ",";
	stainingArray = Array.concat(stainingArray, stainingArea);
	
	
	//***************Morphology Analysis**********************
	//delete pixles outside of ROIs (obtained in extracted MBP image) in original result image
	//selectWindow("Result of " + staining_channelName);
	selectWindow(staining_channelName);
	run("Enhance Contrast...", "saturated=0.5");
	
	//if there was no MBP positive ROI remaining, delete whole image
	if(roiManager("count")==0){
		print("roiManager empty");
		run("Select All");
		setBackgroundColor(0, 0, 0);
		run("Clear", "slice");
	} else{
		//delete pixles outside of ROIs in black and white image
		roiManager("Show All");
		roiManager("Combine");
		run("Clear Outside");
	}

	//calculate morphology of cells only if activated by user
	if(morphologyMeasureBoolean){
		measureMorphologyFunc();
	}
	
	//close not edited original MBP image, so that co-lokalization doesnt get confused. (otherwise two images with same name)
	//close(staining_channelName);
	//rename MBP_black_white to original image name to use it in co-lokalization later
	//selectWindow(staining+"_black_white");
	//rename(staining_channelName);
}

//saves extracted image
function saveExtractedImage(){
	//save copy of image
	run("Duplicate...", " ");
	saveAs("png", outputDirWithOrigImageNameWithStaining +
	"_extracted_image" +".png");
	close();
}

function measureMorphologyFunc(){
		//Select binary black-white image
		//O4
		close(staining + "_black_white");
		selectWindow("Result of " + staining + "_black_white");
		
		//selectWindow(staining_channelName);
		
		//if MBP
		//if(staining == "MBP"){
		//close(staining + "_black_white");
		//selectWindow("Result of " + staining + "_black_white");
		//}else{
		//	selectWindow(staining_channelName);
		//}
	
		//run SKELETONIZE
		skeletonizeImage();
	
		//only add morphology variable names to headings, if needed later used to display result table
		headings += "totalBranchLengthRD" + "," +
		"totalJunctionsRD"+ "," + "totalBranchLengthSkel" + "," + "totalJunctionsSkel"+ "," + "totalTriplePointsSkel" + 
		"," + "totalQuadruplePointsSkel" + ",";
		stainingArray = Array.concat(stainingArray, totalBranchLengthRD, totalJunctionsRD, totalBranchLengthSkel, totalJunctionsSkel, totalTriplePointsSkel, totalQuadruplePointsSkel);
		
		if(staining=="O4"){
			selectWindow("Result of " + staining + "_black_white");
		}else {
			selectWindow(staining_channelName);
		}
		
		//run RIDGE DETECTION
		ridgeDetectionImage(2,1,1.8,0,0.2);
}


function skeletonizeImage(){
	//duplicate image to use it for skeletonize and save black and white image for colokalization later
	run("Duplicate...", " ");
	//run 8bit and Auto threshold since Skeletonize performs error on dublicated image.
	run("8-bit");
	run("Auto Threshold", "method=Huang ignore_black ignore_white white");
	
	//Skeletonize on black-white image
	run("Skeletonize");
	//Analyze Skeleton to obtain measurable results
	run("Analyze Skeleton (2D/3D)", "prune=none calculate show display");
	//rename Result table, including amount of branches, junctions
	Table.rename("Results", staining + "_Skeletonize_Summary");
	//sum up all junctions to total junctions
	for (i=0; i< Table.size; i++){
		totalJunctionsSkel += Table.get("# Junctions", i);
		totalTriplePointsSkel += Table.get("# Triple points", i);
		totalQuadruplePointsSkel += Table.get("# Quadruple points", i);
	}
	////save Table as 
	//saveAs("results", outputDirWithOrigImageNameWithStaining +
	//"_Skeletonize_Summary" +".tsv");
	//close(shortOriginalImageName + "_" + staining +
	//"_Skeletonize_Summary" +".tsv");
	close(staining + "_Skeletonize_Summary");
	
	//rename Branch information table, including amount of branches, junctions
	Table.rename("Branch information", staining + "_Skeletonize_Summary");
	//sum up all branch lengths to total branch length
	for (i=0; i< Table.size; i++){
		totalBranchLengthSkel += Table.get("Branch length", i);
	}
	////save Table as 
	//saveAs("results", outputDirWithOrigImageNameWithStaining +
	//"_Skeletonize_Branch_information" +".tsv");
	//close(shortOriginalImageName + "_" + staining +
	//"_Skeletonize_Branch_information" +".tsv");
	close(staining + "_Skeletonize_Summary");
	
	//close skeletonized image
	//close(shortOriginalImageName+"-3");
	close("Result of "+staining+"_black_white-1");
	//obtained image with longest paths in red of each "island"
	selectWindow("Longest shortest paths");
	close();
	//obtained image with skeleton
	selectWindow("Result-labeled-skeletons");
	close();
	//obtained image
	selectWindow("Tagged skeleton");
	saveAs("png", outputDirWithOrigImageNameWithStaining + 
	"_Skeletonize_Tagged_skeleton" +".png");
	close();
}

function ridgeDetectionImage(highContrast, lowContrast, sigma, lowerTreshold, upperTreshold){
	//close tables before measuring ridge detection
	close("Results"); close("Junctions"); close("Summary");
	close(staining + "_Ridge_Detection_Contour_Class");
	close(staining + "_Ridge_Detection_Junctions");
	close(staining + "_Ridge_Detection_Summary");
	//Ridge Detection for detecting branching morphology of MBP signal
	run("Ridge Detection", "line_width=3.5 high_contrast="+highContrast+ " " +
	"low_contrast="+lowContrast+" displayresults extend_line" +
	"displayresults method_for_overlap_resolution=NONE " +
	"sigma="+sigma+" lower_threshold="+lowerTreshold+" upper_threshold="+upperTreshold+ " "+
	"minimum_line_length=0 maximum=0");
	//save obtained image with skeleton of ridge detection
	saveAs("png", outputDirWithOrigImageNameWithStaining +
	"_Ridge_Detection_obtained_image" +".png");
	close();
	//Create new column in Summary table containing the class of each Contour ID
	//createClassColumnInSummary();
	
	//rename Results table, including coordinates of points which result in contours, length of contour, contours class,
	Table.rename("Results", staining + "_Ridge_Detection_Contour_Class");
	//save Table as 
	//saveAs("results", outputDirWithOrigImageNameWithStaining +
	//"_Ridge_Detection_Contour_Class" +".tsv");
	//close(shortOriginalImageName + "_" + staining +"_Ridge_Detection_Contour_Class" +".tsv");
	close(staining + "_Ridge_Detection_Contour_Class");
	
	//rename Junctions table, including Junction points with X and Y coordinates and Contour IDs
	Table.rename("Junctions", staining + "_Ridge_Detection_Junctions");
	////save Table as 
	//saveAs("results", outputDirWithOrigImageNameWithStaining +
	//"_Ridge_Detection_Junctions" +".tsv");
	//save amount of junctions detected with ridge detection
	totalJunctionsRD = Table.size;
	//close(shortOriginalImageName + "_" + staining +
	//"_Ridge_Detection_Junctions" +".tsv");
	close(staining + "_Ridge_Detection_Junctions");
	
	//rename Summary table, including Contour ID and length of each contour
	Table.rename("Summary", staining + "_Ridge_Detection_Summary");
	////save Table as 
	//saveAs("results", outputDirWithOrigImageNameWithStaining +
	//"_Ridge_Detection_Summary" +".tsv");
	//calculate and save total length of traced paths with ridge
	for (i=0; i< Table.size; i++){
		totalBranchLengthRD += Table.get("Length", i);
	}
	//close(shortOriginalImageName + "_" + staining +
	//"_Ridge_Detection_Summary" +".tsv");
	close(staining + "_Ridge_Detection_Summary");
}


function saveData() {
	directoryWithoutFileSeparator = substring(currentDirectory, 0, lengthOf(currentDirectory)-1);
	lastDirectoryName = substring(directoryWithoutFileSeparator, lastIndexOf(directoryWithoutFileSeparator, "_")+1);
	//Save file in user chosen Output Directory
	fileName = userChosenOutputDirectory + "/" + dataTable + ".tsv";

	headings += "AB_AreaPerCell" + ",";
	if(morphologyMeasureBoolean){
		headings += "branchLengthRDPerCell" + "," + "branchLengthSkelPerCell" + "," + "junctionsRDPerCell" + "," + "junctionsSkelPerCell" + "," + "triplePointsSkelPerCell" + "," + "quadruplePointsSkelPerCell" + ",";
	}
	
	if (endsWith(headings, ",")) {
		headings = substring(headings, 0, lengthOf(headings)-1);
	}
	headingArray = split(headings, ",");
	//referenceObjectCount = Objects = DAPI nuclei 
	// from above: dataArray[0] = referenceObjectCount;
	//originalObjects = Table.getColumn("Objects", "BFE_Results");
	extractedObjectsArray = Table.getColumn("Extracted", "BFE_Results");
	
	//calculate percentage of colocalization of nuclei
	var colocPercentage = newArray(lengthOf(extractedObjectsArray));
	var stainingAreaPerCell = (stainingArea[0]/extractedObjectsArray[0]);
	if(morphologyMeasureBoolean){
		var branchLengthRDPerCell = (totalBranchLengthRD/extractedObjectsArray[0]);
		var branchLengthSkelPerCell = (totalBranchLengthSkel/extractedObjectsArray[0]);
		var junctionsRDPerCell = (totalJunctionsRD/extractedObjectsArray[0]);
		var junctionsSkelPerCell = (totalJunctionsSkel/extractedObjectsArray[0]);
		var triplePointsSkelPerCell = (totalTriplePointsSkel/extractedObjectsArray[0]);
		var quadruplePointsSkelPerCell = (totalQuadruplePointsSkel/extractedObjectsArray[0]);
	}
	
	for (i = 0; i < lengthOf(extractedObjectsArray); i++) {
		colocPercentage[i] = (extractedObjectsArray[i]/referenceObjectCount);
		}
	
	dataArray = Array.concat(shortOriginalImageName, referenceObjectCount, extractedObjectsArray, colocPercentage, stainingArray, stainingAreaPerCell);
	if(morphologyMeasureBoolean){
		meanValueArray = newArray(branchLengthRDPerCell, branchLengthSkelPerCell, junctionsRDPerCell, junctionsSkelPerCell, triplePointsSkelPerCell, quadruplePointsSkelPerCell);
		dataArray = Array.concat(dataArray, meanValueArray);
	}
		
	
	rowCount = Table.size(dataTable);
	for (d = 0; d < headingArray.length; d++) {
		//Table.set(columnName, rowIndex, value, tableName);
		Table.set(headingArray[d], rowCount, dataArray[d], dataTable);
	}
	Table.update(dataTable);
	
	print("fileName"+fileName);
	print("dataTable"+dataTable);
	Table.save(fileName, dataTable);
	close("Results");
	close("Summary");
	close("Log");
	
	//reset headings (delete extra 
	headings = "Image,Nuclei,";
	originalImageName="";
	shortOriginalImageName ="";
	extractedImageName ="";
	outputDirWithOrigImageNameWithStaining="";
	activeChannelCounter = 0;
//var referenceChannelTag = "DAPI";
	referenceChannelName = "";
	referenceChannelID = 0;
	referenceChannelIndex=0;
	referenceObjectCount = 0;
	stainingArray = newArray(0);
	stainingArea= newArray(0);
	totalBranchLengthSkel=0;
	totalBranchLengthRD=0;
	totalJunctionsRD=0;
	totalJunctionsSkel=0;
	totalTriplePointsSkel=0;
	totalQuadruplePointsSkel=0;
}



function getDateStamp() {
	
	getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);

	month += 1;
	if (month < 10) {
		month = "0" + month;
	}
	if (dayOfMonth < 10) {
		dayOfMonth = "0" + dayOfMonth;
	}
	if (hour < 10) {
		hour = "0" + hour;
	}
	if (minute < 10) {
		minute = "0" + minute;
	}
	if (second < 10) {
		second = "0" + second;
	}
	
	dateStamp = ""+ year +""+ month +""+ dayOfMonth +""+ hour +""+ minute +""+ second;
	
	return dateStamp;
}




//currently not functional

/*
function saveDataToExcel() {
	run("Excel Macro Extensions", "debuglogging=false");
	directoryWithoutFileSeparator = substring(currentDirectory, 0, lengthOf(currentDirectory)-1);
	lastDirectoryName = substring(directoryWithoutFileSeparator, lastIndexOf(directoryWithoutFileSeparator, "_")+1);
	fileName = userChosenDirectory + dataTable + ".xlsx";

	if (endsWith(headings, ",")) {
		headings = substring(headings, 0, lengthOf(headings)-1);
	}
	headingArray = split(headings, ",");
	extractedObjectsArray = Table.getColumn("Extracted", "BFE_Results");
	baseArray = newArray(originalImageName, referenceObjectCount);
	dataArray = Array.concat(baseArray, extractedObjectsArray);
	//Array.print(headingArray);
	//Array.print(dataArray);
	
	if (isOpen(dataTable)) {
		for (d = 0; d < headingArray.length; d++) {
			Table.set(headingArray[d], 0, dataArray[d], dataTable);
		}
		Ext.xlsx_SaveTableAsWorksheet(dataTable, fileName, lastDirectoryName, true);
		close(dataTable);
	} else {
		Ext.xlsx_AppendArrayAsExcelRow(dataArray, fileName, lastDirectoryName, 0);
	}
	headings = "Image,Nuclei,";
}

 */