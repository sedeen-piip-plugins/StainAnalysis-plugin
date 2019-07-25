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

#include "algorithm/AlgorithmBase.h"
#include "algorithm/Parameters.h"
#include "algorithm/Results.h"

#include <omp.h>
#include <Windows.h>
#include <fstream>
#include <filesystem> //Requires C++17

//Plugin headers
#include "StainProfile.h"
#include "ColorDeconvolutionKernel.h"

namespace sedeen {
namespace tile {

} // namespace tile

namespace algorithm {
//#define round(x) ( x >= 0.0f ? floor(x + 0.5f) : ceil(x - 0.5f) )

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
	bool buildPipeline(std::shared_ptr<StainProfile>, bool);

    ///Create a text report that combines the output of the stain profile and pixel fraction reports
    std::string generateCompleteReport(std::shared_ptr<StainProfile>) const;
    ///Create a text report summarizing the stain vector profile
	std::string generateStainProfileReport(std::shared_ptr<StainProfile>) const;
    ///Create a text report stating what fraction of the processing area is covered by the filtered output
	std::string generatePixelFractionReport(void) const;

private:
    ///Names of the default stain profile files
    inline static const std::string HematoxylinPEosinSampleFilename()     { return "HematoxylinPEosinSample.xml"; }
    inline static const std::string HematoxylinPEosinFromRJFilename()     { return "HematoxylinPEosinFromRJ.xml"; }
    inline static const std::string HematoxylinPDABFromRJFilename()       { return "HematoxylinPDABFromRJ.xml"; }
    inline static const std::string HematoxylinPEosinPDABFromRJFilename() { return "HematoxylinPEosinPDABFromRJ.xml"; }
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

    OpenFileDialogParameter m_openProfile;
    OptionParameter m_stainSeparationAlgorithm;
    OptionParameter m_stainVectorProfile;
    GraphicItemParameter m_regionToProcess; //single output region

    OptionParameter m_stainToDisplay;
    BoolParameter m_applyThreshold;
    /// User defined Threshold value.
    algorithm::DoubleParameter m_threshold;

    /// The output result
    ImageResult m_result;			
    TextResult m_outputText;
    std::string m_report;

    /// The intermediate image factory after color deconvolution
    std::shared_ptr<image::tile::Factory> m_colorDeconvolution_factory;

    std::ofstream log_file;

private:
    //Member variables
    std::vector<std::string> m_stainVectorProfileOptions;
    std::vector<std::string> m_stainToDisplayOptions;
    double m_thresholdDefaultVal;
    double m_thresholdMaxVal;
};

} // namespace algorithm
} // namespace sedeen

#endif

