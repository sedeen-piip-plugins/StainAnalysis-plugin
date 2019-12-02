/*=============================================================================
 *
 *  Copyright (c) 2019 Sunnybrook Research Institute
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
#include "StainVectorMacenko.h"
#include "StainVectorNiethammer.h"
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
    m_applyThreshold(),
    m_threshold(),
    m_result(),
    m_outputText(),
    m_report(""),
    m_thresholdDefaultVal(20.0),
    m_thresholdMaxVal(300.0),
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

    //Define the default list of names of stains to display
    m_stainToDisplayOptions.push_back("Stain 1");
    m_stainToDisplayOptions.push_back("Stain 2");
    m_stainToDisplayOptions.push_back("Stain 3");
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

	//Get the list of available stain separation algorithms from a temp StainProfile object
	auto tempStainProfile = std::make_shared<StainProfile>();
	std::vector<std::string> tempStainSeparationOptions = tempStainProfile->GetStainSeparationAlgorithmOptions();
	tempStainProfile.reset();
    m_stainSeparationAlgorithm = createOptionParameter(*this, "Stain Separation Algorithm",
        "Select the stain separation algorithm to use to separate the stain components", 0,
		tempStainSeparationOptions, false);

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
    m_applyThreshold = createBoolParameter(*this, "Display with Threshold Applied",
        "If Display with Threshold Applied is set, the threshold value in the slider below will be applied to the stain-separated image",
        true, false); //default value, optional

    // Init the user defined threshold value
	//TEMPORARY: Can't set precision on DoubleParameter right now, so use 1/100 downscale
    auto color = getColorSpace(image);
    auto max_value = (1 << bitsPerChannel(color)) - 1;
    m_threshold = createDoubleParameter(*this,
        "OD x100 Threshold",   // Widget label
        "A Threshold value",   // Widget tooltip
        m_thresholdDefaultVal, // Initial value
        0.0,                   // minimum value
        m_thresholdMaxVal,     // maximum value
        false);

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

		m_result.update(m_colorDeconvolution_factory, m_displayArea, *this);

        //A previous version included an "intermediate result" here, to display
        //a blurry temporary image (rather than just black) while calculations proceeded

		// Update the output text report
		if (false == askedToStop()) {
			auto report = generateCompleteReport(chosenStainProfile);
			
            //TEMPORARY!!!
            //m_outputText.sendText(report);
		}
	}

	// Ensure we run again after an abort
	// a small kludge that causes buildPipeline() to return TRUE
	if (askedToStop()) {
        m_colorDeconvolution_factory.reset();
	}
}//end run

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
         || m_applyThreshold.isChanged() 
         || m_threshold.isChanged() 
         || m_displayArea.isChanged()
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



        //Do this to test out the Macenko method of getting stain vectors
        //Reorganize the UI later


        //Values to use in getting the stain vectors
        auto display_resolution = getDisplayResolution(image(), m_displayArea);

        ////Build the color deconvolution channel
        int numStains = 2;   //TEMPORARY! m_numberOfStainComponents;


        //std::shared_ptr<StainProfile> MacenkoStainProfile = std::make_shared<StainProfile>(*(m_stainProfileList.at(1)));


        if ((numStains <= 0) || (numStains > 3)) {
            return false;
        }
        else if (numStains == 2) {

            //Get the two stain vectors with the Macenko method
            //I also want to know how long this takes to run.


            //here!!! here's the place for trying the next new thing.

            //Create an object to get stain vectors from, using the Macenko algorithm
            //std::shared_ptr<sedeen::image::StainVectorMacenko> stainsFromMacenko 
            //    = std::make_shared<sedeen::image::StainVectorMacenko>(source_factory);
            //stainsFromMacenko->ComputeStainVectors(conv_matrix, 1000, 0.15, 1.0);


            //I need some priors to test with. R+J?
            double priors[9] = { 0.65,0.70,0.29,0.07,0.99,0.11,0.0,0.0,0.0};

            //Create an object to get stain vectors from, using the Niethammer algorithm
            std::shared_ptr<sedeen::image::StainVectorNiethammer> stainsFromNiethammer
                = std::make_shared<sedeen::image::StainVectorNiethammer>(source_factory);
            //stainsFromNiethammer->ComputeStainVectors(conv_matrix, 1000, 0.15, 1.0);


            std::ostringstream ss;

            //Let's try some things.


            //pick up here!!! test the other Set/Get overloads

            stainsFromNiethammer->SetPriors(priors);
            double outPriors[9];

            stainsFromNiethammer->GetPriors(outPriors);

            ss << "So, can I set and get priors treating them as C arrays?" << std::endl;
            for (int i = 0; i < 9; i++) {
                ss << outPriors[i] << ", ";
            }
            ss << std::endl;
          

            //TEMPORARY!
            //std::ostringstream ss;
            //ss << "Here is the output from getting the vectors by Macenko: " << std::endl;
            //ss << "Here is the output from getting the vectors by Niethammer: " << std::endl;
            //for (int i = 0; i < 9; i++) {
            //    ss << conv_matrix[i] << ", ";
            //}


            ss << std::endl;
            m_outputText.sendText(ss.str());


        }
        else {
            m_outputText.sendText("Currently testing algorithms that only accept two stains. Set the number of stains to 2.");
            return false;
        }



        //TEMPORARY
        //Change the contents of the chosenStainProfile
        //REENABLE THIS WHEN YOU'RE READY TO TEST
        //chosenStainProfile->SetProfilesFromDoubleArray(conv_matrix);




        //Scale down the threshold to create more precision
        auto colorDeconvolution_kernel =
            std::make_shared<image::tile::ColorDeconvolution>(DisplayOption, chosenStainProfile, 
                m_applyThreshold, m_threshold/100.0);  //Need to tell it whether to use the threshold or not

        // Create a Factory for the composition of these Kernels
        auto non_cached_factory =
            std::make_shared<FilterFactory>(source_factory, colorDeconvolution_kernel);

        // Wrap resulting Factory in a Cache for speedy results
        m_colorDeconvolution_factory =
            std::make_shared<Cache>(non_cached_factory, RecentCachePolicy(30));

        pipeline_changed = true;
    }//end if parameter values changed

    //
    // Constrain processing to the region of interest provided, if set
    std::shared_ptr<GraphicItemBase> region = m_regionToProcess;
    if (pipeline_changed && (nullptr != region)) {
        // Constrain the output of the pipeline to the region of interest provided
        auto constained_factory = std::make_shared<RegionFactory>(m_colorDeconvolution_factory, region->graphic());

        // Wrap resulting Factory in a Cache for speedy results
        m_colorDeconvolution_factory = std::make_shared<Cache>(constained_factory, RecentCachePolicy(30));
    }

    return pipeline_changed;
}//end buildPipeline

std::string StainAnalysis::generateCompleteReport(std::shared_ptr<StainProfile> theProfile) const {
    //Combine the output of the stain profile report
    //and the pixel fraction report, return the full string
    std::ostringstream ss;
    ss << generateStainProfileReport(theProfile);
    ss << std::endl;
    ss << generatePixelFractionReport();
    return ss.str();
}//end generateCompleteReport

std::string StainAnalysis::generateStainProfileReport(std::shared_ptr<StainProfile> theProfile) const
{
    //If the pointer to theProfile is null, return an error description
    if(theProfile == nullptr) { 
        return "Error reading the stain profile. Please try a different stain profile or restart."; 
    }

    int numStains = theProfile->GetNumberOfStainComponents();
    if (numStains < 0) {
        return "Error reading the stain profile. Please try a different stain profile or restart.";
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
    //Complete, return the string
    return ss.str();
}//end generateStainProfileReport

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
    ss << "stain, above threshold : ";
	ss << std::fixed << std::setprecision(3) << coveredFraction*100  << " %" << std::endl;
	ss << std::endl;

	return ss.str();
}//end generatePixelFractionReport

} // namespace algorithm
} // namespace sedeen
