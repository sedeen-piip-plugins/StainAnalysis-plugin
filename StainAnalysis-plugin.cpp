/*=============================================================================
 *
 *  Copyright (c) 2020 Sunnybrook Research Institute
 *
 *  Permission is hereby granted, free of charge, to any person obtaining a copy
 *  of this software and associated documentation files (the "Software"), to deal
 *  in the Software without restriction, including without limitation the rights
 *  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 *  copies of the Software, and to permit persons to whom the Software is
 *  furnished to do so, subject to the following conditions:
 *
 *  The above copyright notice and this permission notice shall be included in all
 *  copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 *  SOFTWARE.
 *
 *=============================================================================*/

// StainAnalysis-plugin.cpp : Defines the exported functions for the DLL application.
//
#include "StainAnalysis-plugin.h"
#include "ODConversion.h"

#include <sstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

// Sedeen headers
#include "Algorithm.h"
#include "Geometry.h"
#include "Global.h"
#include "Image.h"
#include "image/io/Image.h"
#include "image/tile/Factory.h"

#include <cmrc/cmrc.hpp>

CMRC_DECLARE(stain);

// Poco header needed for the macros below
#include <Poco/ClassLibrary.h>

// Declare that this object has AlgorithmBase subclasses
//  and declare each of those sub-classes
POCO_BEGIN_MANIFEST(sedeen::algorithm::AlgorithmBase)
POCO_EXPORT_CLASS(sedeen::algorithm::StainAnalysis)
POCO_END_MANIFEST

namespace sedeen {
namespace algorithm {

StainAnalysis::StainAnalysis()
	: m_displayArea(),
    m_openProfile(),
    m_stainSeparationAlgorithm(),
    m_stainVectorProfile(),
    m_regionToProcess(),
    m_stainToDisplay(),
    m_applyDisplayThreshold(),
    m_displayThreshold(),
    m_saveSeparatedImage(),
    m_saveFileFormat(),
    m_saveFileAs(),
    m_result(),
    m_outputText(),
    m_report(""),
    m_displayThresholdDefaultVal(20.0),
    m_displayThresholdMaxVal(300.0),
    m_pixelWarningThreshold(1e7), //10,000,000 pixels
    m_colorDeconvolution_factory(nullptr)
{
    // Build the list of stain vector file names
    m_stainProfileFullPathNames.push_back("");  // Leave a blank place for the loaded file
    // HematoxylinPEosinSample
    m_stainProfileFullPathNames.push_back(HematoxylinPEosinSampleFilename());
    // HematoxylinPEosin from Ruifrok and Johnston
    m_stainProfileFullPathNames.push_back(HematoxylinPEosinFromRJFilename());
    // HematoxylinPDAB from Ruifrok and Johnston
    m_stainProfileFullPathNames.push_back(HematoxylinPDABFromRJFilename());
    // HematoxylinPEosinPDAB from Ruifrok and Johnston
    m_stainProfileFullPathNames.push_back(HematoxylinPEosinPDABFromRJFilename());

    //Create stain vector profiles for the loaded and default profiles, push to vector
    //Loaded
    m_LoadedStainProfile = std::make_shared<StainProfile>();
    m_stainProfileList.push_back(m_LoadedStainProfile);
    m_stainVectorProfileOptions.push_back("Loaded From File");
    //Loop over the rest
    for (auto it = m_stainProfileFullPathNames.begin()+1; it != m_stainProfileFullPathNames.end(); ++it) {
        std::shared_ptr<StainProfile> tsp = std::make_shared<StainProfile>();
        //Convert path to string, read stain profile, and check whether the read was successful
        auto const fs = cmrc::stain::get_filesystem();
        if (auto const is_file = fs.is_file((*it).generic_string())) {
            auto const file = fs.open((*it).generic_string());
            if (tsp->readStainProfile(file.begin(), file.size())) {
                m_stainProfileList.push_back(tsp);
                // Get the name of the profile
                m_stainVectorProfileOptions.push_back(tsp->GetNameOfStainProfile());
                continue;
            }
        }

        // make a dummy one
        m_stainProfileList.push_back(std::make_shared<StainProfile>());
        m_stainVectorProfileOptions.push_back("Profile failed to load");
    }

    //Populate the analysis model and separation algorithm lists with a temporary StainProfile
    auto tempStainProfile = std::make_shared<StainProfile>();
    //The stain analysis model options
    m_stainAnalysisModelOptions = tempStainProfile->GetStainAnalysisModelOptions();
    //The stain separation algorithm options
    m_separationAlgorithmOptions = tempStainProfile->GetStainSeparationAlgorithmOptions();
    //Clean up
    tempStainProfile.reset();

    //Define the default list of names of stains to display
    m_stainToDisplayOptions.push_back("Stain 1");
    m_stainToDisplayOptions.push_back("Stain 2");
    m_stainToDisplayOptions.push_back("Stain 3");

    //Choose what format to write the separated images in
    //Define the list of possible save types (flat image vs whole slide image)
    m_saveFileFormatOptions.push_back("Flat image (tif/png/bmp/gif/jpg)");
    //TODO: enable saving a whole slide image
    //m_saveFileFormatOptions.push_back("Whole Slide Image (.svs)");

    //List the actual extensions that should be included in the save dialog window
    m_saveFileExtensionText.push_back("tif");
    m_saveFileExtensionText.push_back("png");
    m_saveFileExtensionText.push_back("bmp");
    m_saveFileExtensionText.push_back("gif");
    m_saveFileExtensionText.push_back("jpg");
    //TODO: enable saving a whole slide image
    //m_saveFileExtensionText.push_back("svs");

}//end constructor

StainAnalysis::~StainAnalysis() {
}//end destructor

void StainAnalysis::init(const image::ImageHandle& image) {
    if (isNull(image)) return;
    //
    // bind algorithm members to UI and initialize their properties

    // Bind system parameter for current view
    m_displayArea = createDisplayAreaParameter(*this);

    //Allow the user to choose a stain vector profile xml file
    sedeen::file::FileDialogOptions openFileDialogOptions = defineOpenFileDialogOptions();
    m_openProfile = createOpenFileDialogParameter(*this, "Stain Profile File",
        "Open a file containing a stain vector profile",
        openFileDialogOptions, true);

    m_stainVectorProfile = createOptionParameter(*this, "Stain Vector Profile",
        "Select the stain vector profile to use; either from the file, or one of the pre-defined profiles", 0,
        m_stainVectorProfileOptions, false);

    m_regionToProcess = createGraphicItemParameter(*this, "Apply to ROI (None for Display Area)",
        "Choose a Region of Interest on which to apply the stain separation algorithm. Choosing no ROI will apply the stain separation to the whole slide image.",
        true); //optional. None means apply to whole slide

    //List of options of the stains to be shown
    m_stainToDisplay = createOptionParameter(*this, "Show Separated Stain",
        "Choose which of the defined stains to show in the display area", 0, m_stainToDisplayOptions, false);

    //User can choose whether to apply the threshold or not
    m_applyDisplayThreshold = createBoolParameter(*this, "Apply Threshold",
        "If Apply Threshold is set, the threshold value in the slider below will be applied to the stain-separated images, including in the saved images",
        true, false); //default value, optional

    // Init the user defined threshold value
	//TEMPORARY: Can't set precision on DoubleParameter right now, so use 1/100 downscale
    //auto color = getColorSpace(image);
    //auto max_value = (1 << bitsPerChannel(color)) - 1;
    m_displayThreshold = createDoubleParameter(*this,
        "OD x100 Threshold",   // Widget label
        "Threshold value to apply to the separated images. Images will be saved with this threshold applied.",   // Widget tooltip
        m_displayThresholdDefaultVal, // Initial value
        0.0,                   // minimum value
        m_displayThresholdMaxVal,     // maximum value
        false);

    //Allow the user to write separated images to file
    m_saveSeparatedImage = createBoolParameter(*this, "Save Separated Image",
        "If checked, the final image will be saved to an output file, of the type chosen in the Save File Format list.",
        false, false);

    m_saveFileFormat = createOptionParameter(*this, "Save File Format",
        "Output image files can be saved as one of five flat image types.",
        0, m_saveFileFormatOptions, false);

    //Allow the user to choose where to save the image files
    sedeen::file::FileDialogOptions saveFileDialogOptions = defineSaveFileDialogOptions();
    m_saveFileAs = createSaveFileDialogParameter(*this, "Save As...",
        "The output image will be saved to this file name. If the file name includes an extension of type TIF/PNG/BMP/GIF/JPG, it will override the Save File Format choice.",
        saveFileDialogOptions, true);

    // Bind result
    m_outputText = createTextResult(*this, "Text Result");
    m_result = createImageResult(*this, " StainAnalysisResult");
      
}//end init

void StainAnalysis::run() {
    //These have to be checked before the parameters are used
    //Check if the stain profile has changed
    bool stainProfile_changed = m_stainVectorProfile.isChanged();
    //Check whether the stain profile file has changed
    bool loadedFile_changed = m_openProfile.isChanged();

    //Get which of the stain vector profiles has been selected by the user
    int chosenProfileNum = m_stainVectorProfile;
    //The user does not have to select a file,
    //but if none is chosen, one of the defaults must be selected
    bool loadResult = LoadStainProfileFromFileDialog();
    //Check whether the user selected to load from file and if so, that it loaded correctly
    if (chosenProfileNum == 0 && (loadResult == false)) {
        m_outputText.sendText("The stain profile file cannot be read. Please click Reset before loading a different file, or choose one of the default profiles.");
        return;
    }
    //else

    // Has display area changed
    bool display_changed = m_displayArea.isChanged();

    //Get the stain profile that should be used
    //Use the vector::at operator to do bounds checking
    std::shared_ptr<StainProfile> chosenStainProfile;
    try {
        chosenStainProfile = m_stainProfileList.at(chosenProfileNum);
    }
    catch (const std::out_of_range& rangeerr) {
        rangeerr.what();
        //The index is out of range. Throw error message
        m_outputText.sendText("The stain profile cannot be found. Choose a default stain profile.");
        return;
    }
    if (chosenStainProfile == nullptr) {
        m_outputText.sendText("The stain profile did not load properly. Click Reset and try another stain profile.");
        return;
    }
    
    //Check the chosenStainProfile to ensure it is built correctly. End running if not
    bool checkProfile = chosenStainProfile->CheckProfile();
    if (!checkProfile) {
        //Something is wrong with the selected stain profile. Show error message, stop running
        m_outputText.sendText("The chosen stain profile did not load properly. Click Reset and try another stain profile.");
        return;
    }

    // Build the operational pipeline
    bool pipeline_changed = buildPipeline(chosenStainProfile, (stainProfile_changed || loadedFile_changed));

	// Update results
	if ( pipeline_changed || display_changed || stainProfile_changed || loadedFile_changed ) {
        //Check whether the user wants to write to image files, that the field is not blank,
        //and that the file can be created or written to
        std::string outputFilePath;
        if (m_saveSeparatedImage == true) {
            //Get the full path file name from the file dialog parameter
            sedeen::algorithm::parameter::SaveFileDialog::DataType fileDialogDataType = this->m_saveFileAs;
            outputFilePath = fileDialogDataType.getFilename();
            //Is the file field blank?
            if (outputFilePath.empty()) {
                m_outputText.sendText("The filename is blank. Please choose a file to save the image to, or uncheck Save Separated Images.");
                return;
            }
            //Does it exist or can it be created, and can it be written to?
            bool validFileCheck = StainProfile::checkFile(outputFilePath, "w");
            if (!validFileCheck) {
                m_outputText.sendText("The file name selected cannot be written to. Please choose another, or check the permissions of the directory.");
                return;
            }
            //Does it have a valid extension? RawImage.save relies on the extension to determine save format
            std::string theExt = getExtension(outputFilePath);
            int extensionIndex = findExtensionIndex(theExt);
            //findExtensionIndex returns -1 if not found
            if (extensionIndex == -1) {
                std::stringstream ss;
                ss << "The extension of the file is not a valid type. The file extension must be: ";
                auto vec = m_saveFileExtensionText;
                for (auto it = vec.begin(); it != vec.end()-1; ++it) {
                    ss << (*it) << ", ";
                }
                std::string last = vec.back();
                ss << "or " << last << ". Choose a correct file type and try again." << std::endl;
                m_outputText.sendText(ss.str());
                return;
            }
        }

        //This is where the magic happens.
        if (nullptr != m_colorDeconvolution_factory) {
            m_result.update(m_colorDeconvolution_factory, m_displayArea, *this);
        }

        //A previous version included an "intermediate result" here, to display
        //a blurry temporary image (rather than just black) while calculations proceeded

        //Determine whether to warn the user about the size of the image to be saved
        bool proceedToSaveImage = false;
        std::string saveInfoAppendToReport = std::string();
        if (m_saveSeparatedImage == true) {
            //Check whether the output image as specified will have more pixels than the given threshold
            bool largeOutputFlag = CheckOutputImageSize(m_pixelWarningThreshold);
            if (largeOutputFlag) {
                //The output file size will be large and will take a long time to save



                //TODO


                //Prompt user for approval to proceed
                bool userApproval = false;





                if (userApproval) {
                    proceedToSaveImage = true;
                }
                else {
                    //User chose not to proceed with saving
                    std::stringstream ss;
                    ss << std::endl <<  "Saving the stain-separated image was cancelled." << std::endl;
                    saveInfoAppendToReport = ss.str();
                }
            }
            else {
                proceedToSaveImage = true;
            }
        }

		// Update the output text report
		if (false == askedToStop()) {
			std::string report = generateCompleteReport(chosenStainProfile);
            if (!saveInfoAppendToReport.empty()) {
                report.append(saveInfoAppendToReport);
            }

            //If an output file should be written and the algorithm ran successfully, save images
            if ((m_saveSeparatedImage == true) && (proceedToSaveImage == true)) {
                //Save the result as a flat image file
                bool saveResult = SaveFlatImageToFile(outputFilePath);
                //Check whether saving was successful
                std::stringstream ss;
                if (saveResult) {
                    ss << std::endl << "Stain-separated image saved as " << outputFilePath << std::endl;
                    report.append(ss.str());
                }
                else {
                    ss << std::endl << "Saving the stain-separated image failed. Please check the file name and directory permissions." << std::endl;
                    report.append(ss.str());
                }
            }

            //Finally, send the report to the results window
            m_outputText.sendText(report);
		}
	}//end if UI changes

	// Ensure we run again after an abort
	// a small kludge that causes buildPipeline() to return TRUE
	if (askedToStop()) {
        m_colorDeconvolution_factory.reset();
	}
}//end run

bool StainAnalysis::buildPipeline(std::shared_ptr<StainProfile> chosenStainProfile, bool somethingChanged) {
    using namespace image::tile;
    bool pipeline_changed = false;

    // Get source image properties
    auto source_factory = image()->getFactory();
    auto source_color = source_factory->getColorSpace();

    //Set this one to change the stain profile
    double conv_matrix[9] = { 0.0 };

    bool doProcessing = false;
    if ( pipeline_changed || somethingChanged
         || m_regionToProcess.isChanged()
         || m_stainSeparationAlgorithm.isChanged() 
         || m_stainVectorProfile.isChanged() 
         || m_stainToDisplay.isChanged() 
         || m_applyDisplayThreshold.isChanged() 
         || m_displayThreshold.isChanged() 
         || m_displayArea.isChanged()
         || m_saveSeparatedImage.isChanged()
         || m_saveFileFormat.isChanged()
         || m_saveFileAs.isChanged()
         || (nullptr == m_colorDeconvolution_factory) )
    {
        //Choose value from the enumeration in ColorDeconvolution
        image::tile::ColorDeconvolution::DisplayOptions DisplayOption;
        switch (m_stainToDisplay)
        {
        case 0:
            DisplayOption = image::tile::ColorDeconvolution::DisplayOptions::STAIN1;
            break;
        case 1:
            DisplayOption = image::tile::ColorDeconvolution::DisplayOptions::STAIN2;
            break;
        case 2:
            DisplayOption = image::tile::ColorDeconvolution::DisplayOptions::STAIN3;
            break;
        default:
            break;
        }

        //Scale down the threshold to create more precision
        auto colorDeconvolution_kernel =
            std::make_shared<image::tile::ColorDeconvolution>(DisplayOption, chosenStainProfile, 
                m_applyDisplayThreshold, m_displayThreshold/100.0);  //Need to tell it whether to use the threshold or not

        // Create a Factory for the composition of these Kernels
        auto non_cached_factory =
            std::make_shared<FilterFactory>(source_factory, colorDeconvolution_kernel);

        // Wrap resulting Factory in a Cache for speedy results
        m_colorDeconvolution_factory =
            std::make_shared<Cache>(non_cached_factory, RecentCachePolicy(30));

        pipeline_changed = true;
    }//end if parameter values changed

    // Constrain processing to the region of interest provided, if set
    std::shared_ptr<GraphicItemBase> region = m_regionToProcess;
    if (pipeline_changed && (nullptr != region)) {
        // Constrain the output of the pipeline to the region of interest provided
        auto constrained_factory = std::make_shared<RegionFactory>(m_colorDeconvolution_factory, region->graphic());

        // Wrap resulting Factory in a Cache for speedy results
        m_colorDeconvolution_factory = std::make_shared<Cache>(constrained_factory, RecentCachePolicy(30));
    }

    return pipeline_changed;
}//end buildPipeline

///Define the open file dialog options outside of init
sedeen::file::FileDialogOptions StainAnalysis::defineOpenFileDialogOptions() {
    sedeen::file::FileDialogOptions theOptions;
    theOptions.caption = "Open stain vector profile: ";
    //theOptions.flags = sedeen::file::FileDialogFlags:: None currently needed
    //theOptions.startDir; //no current preference
    //Define the file type dialog filter
    sedeen::file::FileDialogFilter theDialogFilter;
    theDialogFilter.name = "Stain Vector Profile (*.xml)";
    theDialogFilter.extensions.push_back("xml");
    theOptions.filters.push_back(theDialogFilter);
    return theOptions;
}//end defineOpenFileDialogOptions

///Define the save file dialog options outside of init
sedeen::file::FileDialogOptions StainAnalysis::defineSaveFileDialogOptions() {
    sedeen::file::FileDialogOptions theOptions;
    theOptions.caption = "Save separated images as...";
    //theOptions.flags = sedeen::file::FileDialogFlags:: None currently needed
    //theOptions.startDir; //no current preference
    //Define the file type dialog filter
    sedeen::file::FileDialogFilter theDialogFilter;
    theDialogFilter.name = "Image type";
    //Add extensions in m_saveFileExtensionText to theDialogFilter.extensions 
    //Note: std::copy does not work here
    for (auto it = m_saveFileExtensionText.begin(); it != m_saveFileExtensionText.end(); ++it) {
        theDialogFilter.extensions.push_back(*it);
    }
    theOptions.filters.push_back(theDialogFilter);
    return theOptions;
}//end defineSaveFileDialogOptions

///Access the member file dialog parameter, load into member stain profile
bool StainAnalysis::LoadStainProfileFromFileDialog() {
    //Get the full path file name from the file dialog parameter
    sedeen::algorithm::parameter::OpenFileDialog::DataType fileDialogDataType = this->m_openProfile;
    if (fileDialogDataType.empty()) {
        //There is nothing selected in the file dialog box
        return false;
    }
    //else
    auto profileLocation = fileDialogDataType.at(0);
    std::string theFile = profileLocation.getFilename();

    //Does it exist and can it be read from?
    if (StainProfile::checkFile(theFile, "r")) {
        //Read the stain profile, return false if reading fails        
        bool readFileCheck = m_LoadedStainProfile->readStainProfile(theFile);
        if (readFileCheck) {
            m_stainProfileFullPathNames[0] = theFile;
            return true;
        }
        else {
            m_LoadedStainProfile->ClearProfile();
            m_stainProfileFullPathNames[0] = "";
            return false;
        }
    }
    //else
    return false;
}//end LoadStainProfileFromFileDialog

bool StainAnalysis::CheckOutputImageSize(const double &thresholdVal) {
    bool overThresholdFlag = false;
    //Has a region of interest been set?
    bool roiSet = m_regionToProcess.isUserDefined();
    std::shared_ptr<GraphicItemBase> theRegionOfInterest = m_regionToProcess;
    DisplayRegion displayRegion = m_displayArea;
    auto displayAreaSize = displayRegion.output_size;

    //If a region of interest has been set, use the 
    if (roiSet && theRegionOfInterest != nullptr) {
        Rect rect = containingRect(theRegionOfInterest->graphic());
        double sizeFromRect = static_cast<double>(rect.height())
            * static_cast<double>(rect.width());
        overThresholdFlag = (sizeFromRect > thresholdVal) ? true : false;
    }
    else {
        //No region of interest set. Constrain to display area
        double sizeFromDisplayArea = static_cast<double>(displayAreaSize.height())
            * static_cast<double>(displayAreaSize.width());
        overThresholdFlag = (sizeFromDisplayArea > thresholdVal) ? true : false;
    }

    return overThresholdFlag;
}//end CheckOutputImageSize

bool StainAnalysis::SaveFlatImageToFile(const std::string &p) {
    //It is assumed that error checks have already been performed, and that the type is valid
    //In RawImage::save, the used file format is defined by the file extension.
    //Supported extensions are : .tif, .png, .bmp, .gif, .jpg
    std::string outFilePath = p;
    bool imageSaved = false;
    //Access the output from the output factory
    auto outputFactory = m_colorDeconvolution_factory;
    auto compositor = std::make_unique<image::tile::Compositor>(outputFactory);
    sedeen::image::RawImage outputImage;

    //Has a region of interest been set?
    bool roiSet = m_regionToProcess.isUserDefined();
    std::shared_ptr<GraphicItemBase> theRegion = m_regionToProcess;
    //If a region of interest has been set, constrain output to that area
    if (roiSet && theRegion != nullptr) {
        Rect rect = containingRect(theRegion->graphic());
        //If an ROI is set, output the highest resolution (level 0)
        outputImage = compositor->getImage(0, rect);
    }
    else { 
        //No region of interest set. Constrain to display area
        DisplayRegion region = m_displayArea;
        outputImage = compositor->getImage(region.source_region, region.output_size);
    }
    //Save the outputImage to a file at the given location
    imageSaved = outputImage.save(outFilePath);
    return imageSaved; //true on successful save, false otherwise
}//end SaveFlatImageToFile

const std::string StainAnalysis::getExtension(const std::string &p) {
    namespace fs = std::filesystem; //an alias
    const std::string errorVal = std::string(); //empty
    //Convert the string to a filesystem::path
    fs::path filePath(p);
    //Does the filePath have an extension?
    bool hasExtension = filePath.has_extension();
    if (!hasExtension) { return errorVal; }
    //else
    fs::path ext = filePath.extension();
    return ext.string();
}//end getExtension

const int StainAnalysis::findExtensionIndex(const std::string &x) const {
    const int notFoundVal = -1; //return -1 if not found
    //This method works if the extension has a leading . or not
    std::string theExt(x);
    auto range = std::find(theExt.begin(), theExt.end(), '.');
    theExt.erase(range);
    //Find the extension in the m_saveFileExtensionText vector
    auto vec = m_saveFileExtensionText;
    auto vecIt = std::find(vec.begin(), vec.end(), theExt);
    if (vecIt != vec.end()) {
        ptrdiff_t vecDiff = vecIt - vec.begin();
        int extLoc = static_cast<int>(vecDiff);
        return extLoc;
    }
    else {
        return notFoundVal;
    }
}//end fileExtensionIndex

std::string StainAnalysis::generateCompleteReport(std::shared_ptr<StainProfile> theProfile) const {
    //Combine the output of the stain profile report
    //and the pixel fraction report, return the full string
    std::ostringstream ss;
    ss << generateStainProfileReport(theProfile);
    ss << std::endl;
    //ss << generatePixelFractionReport();
    return ss.str();
}//end generateCompleteReport

std::string StainAnalysis::generateStainProfileReport(std::shared_ptr<StainProfile> theProfile) const
{
    //I think using assert is a little too strong here. Use different error handling.
    assert(nullptr != theProfile);

    int numStains = theProfile->GetNumberOfStainComponents();
    if (numStains < 0) {
        return "Error reading the stain profile. Please change your settings and try again.";
    }
    //Get the profile contents, place in the output stringstream
    std::ostringstream ss;
    ss << std::left << std::setw(5);
    ss << "Using stain profile: " << theProfile->GetNameOfStainProfile() << std::endl;
    ss << "Number of component stains: " << numStains << std::endl;
    ss << std::endl;

    //These are cumulative, not if...else
    //Stain one
    if (numStains >= 1) {
        std::array<double, 3> rgb = theProfile->GetStainOneRGB();
        ss << std::left;
        ss << "Stain 1: " << theProfile->GetNameOfStainOne() << std::endl;
        ss << "R: " << std::setw(10) << std::setprecision(5) << rgb[0] <<
            "G: " << std::setw(10) << std::setprecision(5) << rgb[1] <<
            "B: " << std::setw(10) << std::setprecision(5) << rgb[2] <<
            std::endl;
    }
    //Stain two
    if (numStains >= 2) {
        std::array<double, 3> rgb = theProfile->GetStainTwoRGB();
        ss << std::left;
        ss << "Stain 2: " << theProfile->GetNameOfStainTwo() << std::endl;
        ss << "R: " << std::setw(10) << std::setprecision(5) << rgb[0] <<
            "G: " << std::setw(10) << std::setprecision(5) << rgb[1] <<
            "B: " << std::setw(10) << std::setprecision(5) << rgb[2] <<
            std::endl;
    }
    //Stain three
    if (numStains == 3) {
        std::array<double, 3> rgb = theProfile->GetStainThreeRGB();
        ss << std::left;
        ss << "Stain 3: " << theProfile->GetNameOfStainThree() << std::endl;
        ss << "R: " << std::setw(10) << std::setprecision(5) << rgb[0] <<
            "G: " << std::setw(10) << std::setprecision(5) << rgb[1] <<
            "B: " << std::setw(10) << std::setprecision(5) << rgb[2] <<
            std::endl;
    }
    ss << std::endl;

    //Analysis model and parameters
    std::string analysisModel = theProfile->GetNameOfStainAnalysisModel();
    auto analysisModelParameters = theProfile->GetAllAnalysisModelParameters();
    if (!analysisModel.empty()) {
        ss << "Stain analysis model: " << analysisModel << std::endl;
    }
    if (!analysisModelParameters.empty()) {
        ss << generateParameterMapReport(analysisModelParameters) << std::endl;
    }

    //Separation algorithm and parameters
    std::string separationAlgorithm = theProfile->GetNameOfStainSeparationAlgorithm();
    auto separationAlgorithmParameters = theProfile->GetAllSeparationAlgorithmParameters();
    if (!separationAlgorithm.empty()) {
        ss << "Stain separation algorithm: " << separationAlgorithm << std::endl;
    }
    if (!separationAlgorithmParameters.empty()) {
        ss << generateParameterMapReport(separationAlgorithmParameters) << std::endl;
    }

    //Complete, return the string
    return ss.str();
}//end generateStainProfileReport

std::string StainAnalysis::generateParameterMapReport(std::map<std::string, std::string> p) const {
    std::stringstream ss;
    //Possible parameters are: pTypeNumPixels(), pTypeThreshold(), pTypePercentile(), pTypeHistoBins()
    for (auto it = p.begin(); it != p.end(); ++it) {
        std::string key = it->first;
        std::string val = it->second;
        if (!key.compare(StainProfile::pTypeNumPixels())) {
            ss << "Number of pixels sampled: " << val << std::endl;
        }
        else if (!key.compare(StainProfile::pTypeThreshold())) {
            ss << "Optical Density threshold applied when computing stain vectors: " << val << std::endl;
        }
        else if (!key.compare(StainProfile::pTypePercentile())) {
            ss << "Histogram range percentile: " << val << std::endl;
        }
        else if (!key.compare(StainProfile::pTypeHistoBins())) {
            ss << "Number of histogram bins: " << val << std::endl;
        }
        else {
            //Unknown key, output anyway
            ss << key << ": " << val << std::endl;
        }
    }
    return ss.str();
}//end generateParameterMapReport

std::string StainAnalysis::generatePixelFractionReport() const {
    if (m_colorDeconvolution_factory == nullptr) {
        return "Error accessing the color deconvolution factory. Cannot generate pixel fraction report.";
    }

	using namespace image::tile;

	// Get image from the output factory
	auto compositor = std::make_unique<Compositor>(m_colorDeconvolution_factory);

	DisplayRegion region = m_displayArea;
	auto output_image = compositor->getImage(region.source_region, region.output_size);

	// Get image from the input factory
	auto compositorsource = std::make_unique<Compositor>(image()->getFactory());
	auto input_image = compositorsource->getImage(region.source_region, region.output_size);

	if (m_regionToProcess.isUserDefined()) {
		//myss << m_regionToProcess.isUserDefined() << std::endl;
		std::shared_ptr<GraphicItemBase> roi = m_regionToProcess;
		auto display_resolution = getDisplayResolution(image(), m_displayArea);
		Rect rect = containingRect(roi->graphic());
		output_image = compositor->getImage(rect, region.output_size);
	}

	// Determine number of pixels above threshold
    //Try to count quickly, taking RGB components into account
    unsigned int totalNumPixels = output_image.width()*output_image.height();
    std::vector<unsigned short> pixelSetArray(totalNumPixels,0);
    int iWidth = output_image.width();
    int jHeight = output_image.height();
    int i, j = 0;
	#pragma omp parallel
    {
        #pragma omp for private(i) private(j) //collapse(2) is not available in Visual Studio 2017
        for (i = 0; i < iWidth; ++i) {
            for (j = 0; j < jHeight; ++j) {
                //This relies on implicit conversion from boolean operators to integers 0/1
                pixelSetArray[i*jHeight + j] 
                    = (   output_image.at(i, j, 0).as<uint8_t>()
                       || output_image.at(i, j, 1).as<uint8_t>()
                       || output_image.at(i, j, 2).as<uint8_t>() );
            }
        }
    }
    //Using accumulate means no need for comparison operations to 0/1
    int numPixels = std::accumulate(pixelSetArray.begin(), pixelSetArray.end(), 0);
    double coveredFraction = ((double)numPixels) / ((double)totalNumPixels);

	// Calculate results
	std::ostringstream ss;
	ss << std::left << std::setfill(' ') << std::setw(20);
    //ss << "The absolute number of covered pixels is: " << numPixels << std::endl;
    //ss << "The absolute number of ROI pixels is: " << totalNumPixels << std::endl;
    ss << "Percent of processed region covered by" << std::endl; 
    ss << "stain, above the displayed threshold : ";
	ss << std::fixed << std::setprecision(3) << coveredFraction*100  << " %" << std::endl;
	ss << std::endl;

	return ss.str();
}//end generatePixelFractionReport

} // namespace algorithm
} // namespace sedeen
