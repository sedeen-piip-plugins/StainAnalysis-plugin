/*=========================================================================
 *
 *  Copyright (c) 2019 Sunnybrook Research Institute
 *
 *  License terms pending.
 *
 *=========================================================================*/

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
#include <math.h>

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
    //The new parameters!!!
    m_openProfile(),
    m_stainSeparationAlgorithm(),
    m_stainVectorProfile(),
    m_regionToProcess(),
    m_stainToDisplay(),
    m_applyThreshold(),
    m_threshold()

    //Back to the old parameters!!!
	//m_retainment(),
	//m_displayOptions(),
	//m_result(),
	//m_region_interest(),
	//m_output_text(),
	//m_colorDeconvolution_factory(nullptr)
	//m_threshold_factory(nullptr)
{
    //Define the list of available stain separation algorithms
    m_separationAlgorithmOptions.push_back("Ruifrok Colour Deconvolution");

    m_stainVectorProfileOptions.push_back("Load From File");
    m_stainVectorProfileOptions.push_back("Hematoxylin + Eosin");
    m_stainVectorProfileOptions.push_back("Hematoxylin + DAB");
    m_stainVectorProfileOptions.push_back("Hematoxylin + Eosin + DAB");

    //Define the default list of names of stains to display
    m_stainToDisplayOptions.push_back("Stain 1");
    m_stainToDisplayOptions.push_back("Stain 2");
    m_stainToDisplayOptions.push_back("Stain 3");

}//end constructor

StainAnalysis::~StainAnalysis() {
}//end destructor


void StainAnalysis::run() {
	// Has display area changed
	auto display_changed = m_displayArea.isChanged();

	// Build the segmentation pipeline
	//auto pipeline_changed = buildPipeline();

	// Update results
    /****
	if (pipeline_changed || display_changed ) {

		m_result.update(m_colorDeconvolution_factory, m_displayArea, *this);

		//updateIntermediateResult();

		// Update the output text report
		if (false == askedToStop()) {
			auto report = generateReport();
			m_output_text.sendText(report);
		}
	}

	// Ensure we run again after an abort
	//
	// a small kludge that causes buildPipeline() to return TRUE
	if (askedToStop()) {
		m_colorDeconvolution_factory.reset();
	}
    ****/

}//end run

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


	//Search for the csv file to load the stain matrix
	/*path dir_path = sedeen::image::getSourceDescription(image) + "/sedeen/";
	dir_path.filename();
	dir_path.stem();
	std::string file_name = dir_path.stem().string() + "_StainsFile.csv";
	path path_found = path();*/

	//std::string path_to_image =
	//	image->getMetaData()->get(image::StringTags::SOURCE_DESCRIPTION, 0);
	//const std::string temp_str = path_to_image.substr(path_to_image.find_last_of("/\\") + 1);
	//auto found = path_to_image.find_last_of("/\\") +1;
	//m_path_to_root = path_to_image.substr(0, found);


	//path dir_path = m_path_to_root +"sedeen/";
	//std::string file_name = "StainsFile.csv";
	//path path_found = path();

	//bool thereIsAny = image::serachForfile( dir_path, file_name, path_found );
	//Bind stains selection option list
    /****
	std::vector<std::string> options;
	if(thereIsAny)
	{
		options.push_back("From ROI");
		options.push_back("Hematoxylin + Eosin");
		options.push_back("Hematoxylin + DAB");
		options.push_back("Hematoxylin + Eosin + DAB");
		options.push_back("Load From File");
	}
	else
	{
		options.push_back("From ROI");
		options.push_back("Hematoxylin + Eosin");
		options.push_back("Hematoxylin + DAB");
		options.push_back("Hematoxylin + Eosin + DAB");
	}
	m_retainment = createOptionParameter(*this,
		"Selected Stain",
		"Color-convolution matrix",
		1,        // default selection
		options,
		false); // list of all options
    ****/

	// Bind display options for stain components
    /****
    std::vector<std::string> displays(3);
	displays[0] = "STAIN 1";
	displays[1] = "STAIN 2";
	displays[2] = "STAIN 3";
	m_displayOptions = createOptionParameter(*this,
		"Display",
		"Display stain components",
		0,        // default selection
		displays,
		false); // list of all options
    ****/


    // Bind result
	m_output_text = createTextResult(*this, "Text Result");
	m_result = createImageResult(*this, " StainAnalysisResult");

}//end init

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
			//	m_path_to_root = path_to_image.substr(0, found);
			//	m_path_to_stainfile = openFile(m_path_to_root);
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

			image::tile::ColorDeconvolution::DisplyOptions displyOption;
			switch (m_displayOptions)
			{
			case 0:
				displyOption = image::tile::ColorDeconvolution::DisplyOptions::STAIN1;
				break;
			case 1:
				displyOption = image::tile::ColorDeconvolution::DisplyOptions::STAIN2;
				break;
			case 2:
				displyOption = image::tile::ColorDeconvolution::DisplyOptions::STAIN3;
				break;
			default:
				break;
			}

				auto colorDeconvolution_kernel =
					std::make_shared<image::tile::ColorDeconvolution>(retainment, displyOption, conv_matrix, m_threshold, m_path_to_root);

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

std::string StainAnalysis::generateReport(double conv_matrix[9]) const
{
	assert(nullptr != conv_matrix);

	// Format output with 10 character wide columns
	std::ostringstream ss;
	ss << std::left << std::setw(5);
	ss << "Color deconvolution - Stains Coefficients " << std::endl;
	ss << std::endl;

	// Calculate results
	for(int i =0 ;i < 3; ++i) {

		ss << std::left;
		ss << "Stain-" << i+1 << std::endl <<
			"R" << i+1 << ":" << std::setw(10) << std::setprecision(5) << conv_matrix[i*3] <<
			"G" << i+1 << ":" << std::setw(10) << std::setprecision(5) << conv_matrix[i*3+1] <<
			"B" << i+1 << ":" << std::setw(10) << std::setprecision(5) << conv_matrix[i*3+2] ;
		ss << std::endl;
	}

	return ss.str();
}


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






std::string StainAnalysis::openFile(std::string path)
{
    OPENFILENAME ofn;
    char szFileName[MAX_PATH] = "";

    //HACK
    std::string m_path_to_root = "";

    ZeroMemory(&ofn, sizeof(ofn));
    ofn.lStructSize = sizeof(ofn);
    ofn.hwndOwner = NULL;
    ofn.lpstrFilter = "*.csv";
    ofn.lpstrFile = (LPSTR)szFileName;
    ofn.nMaxFile = MAX_PATH;
    ofn.Flags = OFN_EXPLORER | OFN_FILEMUSTEXIST | OFN_HIDEREADONLY;
    //ofn.lpstrDefExt = (LPSTR)L"tif";
    ofn.lpstrInitialDir = (LPSTR)m_path_to_root.c_str();
    GetOpenFileName(&ofn);

    return ofn.lpstrFile;
}//end openFile

} // namespace algorithm
} // namespace sedeen


