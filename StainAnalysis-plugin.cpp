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

// StainAnalysis.cpp : Defines the exported functions for the DLL application.
//
// Primary header
#include "StainAnalysis-plugin.h"

#include <algorithm>
#include <cassert>
#include <sstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>

// DPTK headers
#include "Algorithm.h"
#include "Geometry.h"
#include "Global.h"
#include "Image.h"
#include "image/io/Image.h"
#include "image/tile/Factory.h"

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
    m_pathToPlugin(""),

    //Back to the old parameters!!!
	//m_retainment(),
	//m_displayOptions(),
	m_result(),
	//m_region_interest(),
	m_outputText(),
	m_colorDeconvolution_factory(nullptr),
	m_threshold_factory(nullptr)
{
    //Get the path of the current plugin
    sedeen::file::Location pluginLocation = sedeen::file::getWorkingDirectory();
    if (GetPluginRelativeDirectory() != "PLUGIN_RELATIVE_DIR-NOTFOUND") {
        //Convert to std::filesystem::path type, and use the path concatenation operator
        m_pathToPlugin = std::filesystem::path(pluginLocation.getFilename()) / GetPluginRelativeDirectory();
        //Build the list of stain vector file names
        m_stainProfileFullPathNames.push_back(""); //Leave a blank place for the loaded file
        //HematoxylinPEosin
        m_stainProfileFullPathNames.push_back(m_pathToPlugin / HematoxylinPEosinFilename());
        //HematoxylinPDAB
        m_stainProfileFullPathNames.push_back(m_pathToPlugin / HematoxylinPDABFilename());
        //HematoxylinPEosinPDAB
        m_stainProfileFullPathNames.push_back(m_pathToPlugin / HematoxylinPEosinPDABFilename());        
    }
    else { //Could not get the plugin's relative directory
        //Location of the default stain profiles is unknown
        m_pathToPlugin = "";
        //Do not try to load the default stain profiles
        m_stainProfileFullPathNames.push_back("");
    }

    //Create stain vector profiles for the loaded and default profiles, push to vector
    //Loaded
    m_LoadedStainProfile = std::make_shared<StainProfile>();
    m_stainProfileList.push_back(m_LoadedStainProfile);
    m_stainVectorProfileOptions.push_back("Loaded From File");
    //Loop over the rest
    for (auto it = m_stainProfileFullPathNames.begin()+1; it != m_stainProfileFullPathNames.end(); ++it) {
        std::shared_ptr<StainProfile> tsp = std::make_shared<StainProfile>();
        //Convert path to string, read stain profile, and check whether the read was successful
        bool success = tsp->readStainProfile((*it).generic_string());
        if (success) {
            m_stainProfileList.push_back(tsp);
            //Get the name of the profile
            m_stainVectorProfileOptions.push_back(tsp->GetNameOfStainProfile());
        }
        else {
            //make a dummy one
            m_stainProfileList.push_back(std::make_shared<StainProfile>());
            m_stainVectorProfileOptions.push_back("Profile failed to load");
        }
        tsp.reset();
    }

    //Define the list of available stain separation algorithms
    m_separationAlgorithmOptions.push_back("Ruifrok+Johnston Deconvolution");

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
        "Open a file containing a stain vector profile.",
        openFileDialogOptions, true);

    m_stainSeparationAlgorithm = createOptionParameter(*this, "Stain Separation Algorithm",
        "Select the stain separation algorithm to use to separate the stain components", 0,
        m_separationAlgorithmOptions, false);

    m_stainVectorProfile = createOptionParameter(*this, "Stain Vector Profile",
        "Select the stain vector profile to use; either from the file, or one of the pre-defined profiles", 0,
        m_stainVectorProfileOptions, false);

    m_regionToProcess = createGraphicItemParameter(*this, "Apply to ROI (None for Whole Slide)",
        "Choose a Region of Interest on which to apply the stain separation algorithm. Choosing no ROI will apply the stain separation to the whole slide image.",
        true); //optional. None means apply to whole slide

    //List of options of the stains to be shown
    m_stainToDisplay = createOptionParameter(*this, "Show Separated Stain",
        "Choose which of the defined stains to show in the display area", 0, m_stainToDisplayOptions, false);

    //User can choose whether to apply the threshold or not
    m_applyThreshold = createBoolParameter(*this, "Display with Threshold Applied",
        "If Display with Threshold Applied is set, the threshold value in the slider below will be applied to the stain-separated image",
        false, false); //default value, optional

    // Init the user defined threshold value
    auto color = getColorSpace(image);
    auto max_value = (1 << bitsPerChannel(color)) - 1;
    m_threshold = createDoubleParameter(*this,
        "Threshold", // Widget label
        "A Threshold value", // Widget tooltip
        1.0,         // Initial value
        0.0,          // minimum value
        50.0,        // maximum value
        false);

    //std::string pathToImage =
    //    image->getMetaData()->get(image::StringTags::SOURCE_DESCRIPTION, 0);
    //const std::string temp_str = pathToImage.substr(pathToImage.find_last_of("/\\") + 1);
    //auto found = pathToImage.find_last_of("/\\") + 1;
    //m_pathToRoot = pathToImage.substr(0, found);

    //std::filesystem::path dirPath = m_pathToRoot / "sedeen/";
    //std::filesystem::path pathFound = std::filesystem::path();

    // Bind result
    m_outputText = createTextResult(*this, "Text Result");
    m_result = createImageResult(*this, " StainAnalysisResult");

}//end init


void StainAnalysis::run() {
    //The user does not have to select a file,
    //but if none is chosen, one of the defaults must be selected
    bool loadResult = LoadStainProfileFromFileDialog();
    //Get which of the stain vector profiles has been selected by the user
    int chosenProfileNum = m_stainVectorProfile;

    //Define a pointer for the stain profile to be 
    //Check whether the user selected to load from file and if so, that it loaded correctly
    if (!loadResult && chosenProfileNum == 0) {
        throw std::runtime_error("The stain profile file cannot be read. Choose a default stain profile or try loading a different file.");
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
        throw std::runtime_error("The stain profile cannot be found. Choose a default stain profile.");
        return;
    }
    
    // Build the operational pipeline
    bool pipeline_changed = buildPipeline(chosenStainProfile);

	// Update results
	if ( pipeline_changed || display_changed ) {

		m_result.update(m_colorDeconvolution_factory, m_displayArea, *this);
        //*********************************************************************************************************************************here
		//updateIntermediateResult();

		// Update the output text report
		if (false == askedToStop()) {
			auto report = generateStainAnalysisReport(chosenStainProfile);
			m_outputText.sendText(report);
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
        m_stainProfileFullPathNames[0] = theFile;
        //Clear the loaded stain profile
        //something here

        m_LoadedStainProfile->readStainProfile(theFile);
        return true;
    }
    else {
        return false;
    }
}//end LoadStainProfileFromFileDialog


bool StainAnalysis::buildPipeline(std::shared_ptr<StainProfile> chosenStainProfile) {
    using namespace image::tile;
    bool pipeline_changed = false;
    //There are two possible behaviors: RegionOfInterest, and LoadFromProfile
    //This plugin only uses LoadFromProfile
    image::tile::ColorDeconvolution::Behavior behaviorType = image::tile::ColorDeconvolution::Behavior::LoadFromProfile;

    // Get source image properties
    auto source_factory = image()->getFactory();
    auto source_color = source_factory->getColorSpace();

    bool doProcessing = false;
    //bool regionsChanged(false);
    //Have any of the regions to be processed changed?
    //for (auto it = m_regionsToProcess.begin(); it != m_regionsToProcess.end(); ++it) {
    //    regionsChanged = regionsChanged || (*it).isChanged();
    //}
    //Figure out whether any settings have changed
    //|| regionsChanged
    if ( pipeline_changed 
         || m_regionToProcess.isChanged()
         || m_stainSeparationAlgorithm.isChanged() 
         || m_stainVectorProfile.isChanged() 
         || m_stainToDisplay.isChanged() 
         || m_applyThreshold.isChanged() 
         || m_threshold.isChanged() 
         || m_displayArea.isChanged()
         || (nullptr == m_colorDeconvolution_factory) )
    {
        //Build the color deconvolution channel

        if (behaviorType == image::tile::ColorDeconvolution::Behavior::RegionOfInterest) {
            //for now, do nothing here.

            /**
            if (retainment == image::tile::ColorDeconvolution::Behavior::RegionOfInterest) {

                std::vector<std::shared_ptr<GraphicItemBase>> region_of_interests;

                if(m_region_interest.at(0).isUserDefined() && m_region_interest.at(1).isUserDefined()
                    && m_region_interest.at(2).isUserDefined() )
                {
                    auto display_resolution = getDisplayResolution(image(), m_displayArea);
                    for(int i=0; i < 3; i++)
                    {
                        std::shared_ptr<GraphicItemBase> region = m_region_interest.at(i);
                        region_of_interests.push_back( region ); //region->graphic()

                        Rect rect = containingRect(region_of_interests.at(i)->graphic());
                    }

                    image::getStainsComponents(source_factory,
                        region_of_interests,
                        display_resolution, conv_matrix);

                    m_report = generateReport(conv_matrix);
                }
                else
                {
                    doProcessing = false;
                    retainment = image::tile::ColorDeconvolution::Behavior::HematoxylinPEosin;

                }
            }
        **/

        }//end Behavior is RegionOfInterest

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

        auto colorDeconvolution_kernel =
            std::make_shared<image::tile::ColorDeconvolution>(behaviorType, DisplayOption, chosenStainProfile, m_applyThreshold, m_threshold);  //Need to tell it whether to use the threshold or not

        // Create a Factory for the composition of these Kernels
        auto non_cached_factory =
            std::make_shared<FilterFactory>(source_factory, colorDeconvolution_kernel);

        // Wrap resulting Factory in a Cache for speedy results
        m_colorDeconvolution_factory =
            std::make_shared<Cache>(non_cached_factory, RecentCachePolicy(30));

        pipeline_changed = true;

    }//end if parameter values changed

    //
    // Constrain processing to the region of interest provided
    std::shared_ptr<GraphicItemBase> region = m_regionToProcess;
    if (pipeline_changed && (nullptr != region)) {
        // Constrain the output of the pipeline to the region of interest provided
        auto constained_factory = std::make_shared<RegionFactory>(m_colorDeconvolution_factory, region->graphic());

        // Wrap resulting Factory in a Cache for speedy results
        m_colorDeconvolution_factory = std::make_shared<Cache>(constained_factory, RecentCachePolicy(30));
    }

    return pipeline_changed;
}//end buildPipeline






/****
void StainAnalysis::updateIntermediateResult()
{
	// Update UI with the results of the given factory
	auto update_result = [&](const std::shared_ptr<image::tile::Factory> &factory) {
		// Create a compositor
		auto compositor = std::unique_ptr<image::tile::Compositor>(new image::tile::Compositor(factory));

		// Extract image from it
		//auto source_region = image()->getFactory()->getLevelRegion(0);
		//auto image = compositor->getImage(source_region, Size(source_region.width(), source_region.height()));

		DisplayRegion region = m_displayArea;
		auto image = compositor->getImage(region.source_region, region.output_size);

		// Update UI
		m_result.update(image, region.source_region);
	};

	update_result(m_colorDeconvolution_factory);
}
****/

std::string StainAnalysis::generateStainAnalysisReport(std::shared_ptr<StainProfile> theProfile) const
{
    //I think using assert is a little too strong here. Use different error handling.
	assert(nullptr != theProfile);

    int numStains = theProfile->GetNumberOfStainComponents();
    if (numStains < 0) {
        throw std::runtime_error("Error reading the stain profile. Please try a different stain profile or restart.");
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
}//end generateStainAnalysisReport


std::string StainAnalysis::generateReport() const{

	assert(nullptr != m_colorDeconvolution_factory);

	using namespace image::tile;

	// Get image from the output factory
	auto compositor =
		std::unique_ptr<Compositor>(new Compositor(m_colorDeconvolution_factory));

	DisplayRegion region = m_displayArea;
	auto output_image = compositor->getImage(region.source_region, region.output_size);

	// Get image from the input factory
	auto compositorsource =
		std::unique_ptr<Compositor>(new Compositor(image()->getFactory()));
	auto input_image = compositorsource->getImage(region.source_region, region.output_size);

	if (m_regionToProcess.isUserDefined()) {
		//myss << m_regionToProcess.isUserDefined() << std::endl;
		std::shared_ptr<GraphicItemBase> roi = m_regionToProcess;
		auto display_resolution = getDisplayResolution(image(), m_displayArea);
		Rect rect = containingRect(roi->graphic());
		output_image = compositor->getImage(rect, region.output_size);
	}

	// Determine number of pixels above threshold
	unsigned int totalnumPixels = output_image.width()*output_image.height();
	unsigned int numPixels=0;
	int j = 0;
	#pragma omp parallel for private(j)
	for (int i=0;i<output_image.width();i++){
		for (int j=0;j<output_image.height();j++){
			if(!(output_image.at(i,j, 0).as<uint8_t>() == 0))
			{
				numPixels++;
			}
		}
	}

	// Get pixel size of a pixel in the native resolution.
	Size full_size = getDimensions(image(), 0);
	float scalew = 1.0;
	float scaleh = 1.0;
	unsigned int _totalnumPixels = 1;
	float _numPixels= (float)numPixels;
	if(full_size.width() > output_image.width()){
		scalew = (float)full_size.width() / (float)output_image.width();
		scaleh = (float)full_size.height() / (float)output_image.height();
		_totalnumPixels = full_size.width()*full_size.height();
		_numPixels = (float)(_numPixels*scalew*scaleh)/(float)_totalnumPixels;
		//_totalnumPixels = full_size.width()*full_size.height() - _numPixels;
	}
	else{
		scalew = float(output_image.width())/float(full_size.width());
		scaleh = (float)output_image.height()/(float)full_size.height();

		_totalnumPixels = full_size.width()*full_size.height();
		_numPixels = (float)numPixels/(scalew*scaleh);
		_numPixels = _numPixels/(float)_totalnumPixels;
		//_totalnumPixels = full_size.width()*full_size.height() - _numPixels;
	}

	// Calculate results
	std::ostringstream ss;
	ss << std::left << std::setfill(' ') << std::setw(20);
	ss << "Pixels belong to FG:";
	ss << std::fixed << std::setprecision(3) << _numPixels*100  << " \%" << std::endl;
	ss << std::endl;

	return ss.str();
}



} // namespace algorithm
} // namespace sedeen


/****
bool StainAnalysis::buildPipeline() {
    using namespace image::tile;
    bool pipeline_changed = false;

    // Get source image properties
    auto source_factory = image()->getFactory();
    auto source_color = source_factory->getColorSpace();
    double conv_matrix[9]= {0.0};
    bool ROIIsdefined = (m_region_interest.at(0).isUserDefined() &&
        m_region_interest.at(1).isUserDefined() &&
        m_region_interest.at(2).isUserDefined() ) &&
        (m_region_interest.at(0).isChanged() ||
        m_region_interest.at(1).isChanged() ||
        m_region_interest.at(2).isChanged() );

    bool doProcessing = false;
    if (m_retainment.isChanged() || m_displayOptions.isChanged() ||
        ROIIsdefined || m_threshold.isChanged() || m_regionToProcess.isChanged()  ||
        m_displayArea.isChanged() ||
        (nullptr == m_colorDeconvolution_factory) || pipeline_changed) {
            // Build Color Deconvolution Kernel
            image::tile::ColorDeconvolution::Behavior retainment;
            switch (m_retainment)
            {
            case 0:
                retainment = image::tile::ColorDeconvolution::Behavior::RegionOfInterest;
                break;
            case 1:
                retainment = image::tile::ColorDeconvolution::Behavior::HematoxylinPEosin;
                break;
            case 2:
                retainment = image::tile::ColorDeconvolution::Behavior::HematoxylinPDAB;
                break;
            case 3:
                retainment = image::tile::ColorDeconvolution::Behavior::HematoxylinPEosinPDAB;
                break;
            case 4:
                retainment = image::tile::ColorDeconvolution::Behavior::LoadFromFile;
                break;
            default:
                break;
            }

            //if(retainment == image::tile::ColorDeconvolution::Behavior::LoadFromFile)
            //{
            //	std::string path_to_image =
            //		image()->getMetaData()->get(image::StringTags::SOURCE_DESCRIPTION, 0);
            //	auto found = path_to_image.find_last_of(".");
            //	m_pathToRoot = path_to_image.substr(0, found);
            //	m_path_to_stainfile = openFile(m_pathToRoot);
            //}

            if (retainment == image::tile::ColorDeconvolution::Behavior::RegionOfInterest) {

                std::vector<std::shared_ptr<GraphicItemBase>> region_of_interests;

                if(m_region_interest.at(0).isUserDefined() && m_region_interest.at(1).isUserDefined()
                    && m_region_interest.at(2).isUserDefined() )
                {
                    auto display_resolution = getDisplayResolution(image(), m_displayArea);
                    for(int i=0; i < 3; i++)
                    {
                        std::shared_ptr<GraphicItemBase> region = m_region_interest.at(i);
                        region_of_interests.push_back( region ); //region->graphic()

                        Rect rect = containingRect(region_of_interests.at(i)->graphic());
                    }

                    image::getStainsComponents(source_factory,
                        region_of_interests,
                        display_resolution, conv_matrix);

                    m_report = generateReport(conv_matrix);
                }
                else
                {
                    doProcessing = false;
                    retainment = image::tile::ColorDeconvolution::Behavior::HematoxylinPEosin;

                }
            }

            //m_report = generateReport(conv_matrix);

            image::tile::ColorDeconvolution::DisplayOptions DisplayOption;
            switch (m_displayOptions)
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

                auto colorDeconvolution_kernel =
                    std::make_shared<image::tile::ColorDeconvolution>(retainment, DisplayOption, conv_matrix, m_threshold, m_pathToRoot);

                // Create a Factory for the composition of these Kernels
                auto non_cached_factory =
                    std::make_shared<FilterFactory>(source_factory,  colorDeconvolution_kernel);

                // Wrap resulting Factory in a Cache for speedy results
                m_colorDeconvolution_factory =
                    std::make_shared<Cache>(non_cached_factory, RecentCachePolicy(30));

                pipeline_changed = true;
    }

    //
    // Constrain processing to the region of interest provided
    //
    //
    std::shared_ptr<GraphicItemBase> region = m_regionToProcess;
    if (pipeline_changed && (nullptr != region)) {
        // Constrain the output of the pipeline to the region of interest provided
        auto constained_factory =
            std::make_shared<RegionFactory>(m_colorDeconvolution_factory,
            region->graphic());

        // Wrap resulting Factory in a Cache for speedy results
        m_colorDeconvolution_factory =
            std::make_shared<Cache>(constained_factory, RecentCachePolicy(30));
    }

    return pipeline_changed;
}
****/