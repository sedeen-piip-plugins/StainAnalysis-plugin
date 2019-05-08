/*=========================================================================
 *
 *  Copyright (c) 2019 Sunnybrook Research Institute
 *
 *  License terms pending.
 *
 *=========================================================================*/

#ifndef SEDEEN_SRC_PLUGINS_STAINANALYSIS_STAINANALYSIS_H
#define SEDEEN_SRC_PLUGINS_STAINANALYSIS_STAINANALYSIS_H

//Convert a macro defined in CMake to a string
#define StringLiteral(stl) #stl
#define MacroToString(macro) StringLiteral(macro)
//Specifically, convert PLUGIN_RELATIVE_DIR to a string literal
#ifndef PLUGIN_RELATIVE_DIR
#define PLUGIN_RELATIVE_DIR_STRING "PLUGIN_RELATIVE_DIR-NOTFOUND"
#else
#define PLUGIN_RELATIVE_DIR_STRING MacroToString(PLUGIN_RELATIVE_DIR)
#endif

// DPTK headers - a minimal set 
#include "algorithm/AlgorithmBase.h"
#include "algorithm/Parameters.h"
#include "algorithm/Results.h"

#include "ColorDeconvolutionKernel.h"

#include <omp.h>
#include <Windows.h>
#include <fstream>
#include <filesystem> //Requires C++17

//Plugin headers
#include "StainProfile.h"

namespace sedeen {
namespace tile {

} // namespace tile

namespace algorithm {
#define round(x) ( x >= 0.0f ? floor(x + 0.5f) : ceil(x - 0.5f) )

/// Stain Analysis
/// This plugin implements stain separation using the colour deconvolution
/// method described in:
//Ruifrok AC, Johnston DA. Quantification of histochemical
//staining by color deconvolution. Analytical & Quantitative
//Cytology & Histology 2001; 23: 291-299.
class StainAnalysis : public algorithm::AlgorithmBase {
public:
	StainAnalysis();
    virtual ~StainAnalysis();

private:
	// virtual function
	virtual void run();
	virtual void init(const image::ImageHandle& image);

    //Define the open file dialog options outside of init
    sedeen::file::FileDialogOptions defineOpenFileDialogOptions();

	/// Creates the Color Deconvolution pipeline with a cache
	//
	/// \return 
	/// TRUE if the pipeline has changed since the call to this function, FALSE
	/// otherwise
	bool buildPipeline(std::shared_ptr<StainProfile>);

    ///Create a report including the stain profile data used and the resulting stain analysis output
	std::string generateStainAnalysisReport(std::shared_ptr<StainProfile>) const;
    ///Comment pending
	std::string generateReport(void) const;

	//void StainAnalysis::updateIntermediateResult();

private:
    ///Names of the default stain profile files
    inline static const std::string HematoxylinPEosinFilename()     { return "HematoxylinPEosin.xml"; }
    inline static const std::string HematoxylinPDABFilename()       { return "HematoxylinPDAB.xml"; }
    inline static const std::string HematoxylinPEosinPDABFilename() { return "HematoxylinPEosinPDAB.xml"; }
    ///Returns the directory of the plugin relative to the Sedeen Viewer directory as a string
    inline static const std::filesystem::path GetPluginRelativeDirectory() { return std::filesystem::path(PLUGIN_RELATIVE_DIR_STRING); }

    ///Access the member file dialog parameter, if possible load the stain profile, return true on success
    bool LoadStainProfileFromFileDialog();

    ///std::filesystem::path type to the plugin's directory
    std::filesystem::path m_pathToPlugin;

    ///List of the full path file names of the stain profiles
    std::vector<std::filesystem::path> m_stainProfileFullPathNames;

    ///List of the connected stain profile objects
    std::vector<std::shared_ptr<StainProfile>> m_stainProfileList;
    ///Keep a pointer directly to the loaded stain profile
    std::shared_ptr<StainProfile> m_LoadedStainProfile;

private:
	DisplayAreaParameter m_displayArea;

    //The new parameters!!!
    OpenFileDialogParameter m_openProfile;
    OptionParameter m_stainSeparationAlgorithm;
    OptionParameter m_stainVectorProfile;
    GraphicItemParameter m_regionToProcess; //ONE output region for now
    //std::vector<GraphicItemParameter> m_regionsToProcess;

    OptionParameter m_stainToDisplay;
    BoolParameter m_applyThreshold;
    /// User defined Threshold value.
    algorithm::DoubleParameter m_threshold;
    //End of the new parameters

    /// The output result
    ImageResult m_result;			
    TextResult m_outputText;
    std::string m_report;

    /// The intermediate image factory after color deconvolution
    std::shared_ptr<image::tile::Factory> m_colorDeconvolution_factory;

    /// The image factory after thresholding
    std::shared_ptr<image::tile::Factory> m_threshold_factory;
    std::ofstream log_file;

    //Old parameters
	/// The output result
	/// Parameter for selecting which of the intermediate result to display
	//algorithm::OptionParameter m_output_option;

    /// User region of interest (these are for the ROIs defining reference stains. Now only needed in CreateStainVectorProfile)
	//std::vector<algorithm::GraphicItemParameter> m_region_interest;
	//


private:
    //Member variables
    std::vector<std::string> m_separationAlgorithmOptions;
    std::vector<std::string> m_stainVectorProfileOptions;
    std::vector<std::string> m_stainToDisplayOptions;

};

} // namespace algorithm
} // namespace sedeen

#endif

