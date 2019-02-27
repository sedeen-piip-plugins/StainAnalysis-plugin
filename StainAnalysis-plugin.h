/*=========================================================================
 *
 *  Copyright (c) 2019 Sunnybrook Research Institute
 *
 *  License terms pending.
 *
 *=========================================================================*/

#ifndef SEDEEN_SRC_PLUGINS_STAINANALYSIS_STAINANALYSIS_H
#define SEDEEN_SRC_PLUGINS_STAINANALYSIS_STAINANALYSIS_H

// DPTK headers - a minimal set 
#include "algorithm/AlgorithmBase.h"
#include "algorithm/Parameters.h"
#include "algorithm/Results.h"
#include "ColorDeconvolutionKernel.h"

#include <omp.h>
#include <Windows.h>
#include <fstream>

namespace sedeen {
namespace tile {

} // namespace tile

namespace algorithm {
/// Color Deconvolution
/// This plugin implements stain separation using the colour deconvolution
/// method described in:
//Ruifrok AC, Johnston DA. Quantification of histochemical
//staining by color deconvolution. Analytical & Quantitative
//Cytology & Histology 2001; 23: 291-299.
#define round(x) ( x >= 0.0f ? floor(x + 0.5f) : ceil(x - 0.5f) )

class StainAnalysis : public algorithm::AlgorithmBase {
public:
	StainAnalysis();

private:
	// virtual function
	virtual void run();
	virtual void init(const image::ImageHandle& image);

	/// Creates the Color Deconvolution pipeline with a cache
	//
	/// \return 
	/// TRUE if the pipeline has changed since the call to this function, FALSE
	/// otherwise
	bool buildPipeline();

	/// Generates a report of the stain matrix calculated from ROI
	//
	/// Report is formatted as a table containing the cofficients
	/// of each stain
	//
	/// \return
	/// a string containing the cofficients of each stain calculated from ROI
	std::string generateReport(double[9]) const;
	std::string generateReport(void) const;
	std::string openFile(std::string path);


	void StainAnalysis::updateIntermediateResult();


	/// Search folder to find a specific file containg saved stain matrix
	//
	/// Report is formatted as a table containing the cofficients
	/// of each stain
	//
	/// \return
	/// a string containing the cofficients of each stain calculated from ROI
	//bool find_file( const boost::filesystem::path & , const std::string & , boost::filesystem::path & );        

private:
	std::string m_path_to_root;
	std::string m_path_to_stainfile;

	algorithm::DisplayAreaParameter m_display_area;
	algorithm::OptionParameter m_retainment;
	algorithm::OptionParameter m_displayOptions;
	/// Parameter for selecting threshold retainment 
	//algorithm::OptionParameter m_behavior;
	/// User defined Threshold value.
	algorithm::DoubleParameter m_threshold;
	/// The output result
	algorithm::ImageResult m_result;			
	algorithm::TextResult m_output_text;
	std::string m_report;
	/// Parameter for selecting which of the intermediate result to display
	//algorithm::OptionParameter m_output_option;
	/// User region of interest
	std::vector<algorithm::GraphicItemParameter> m_region_interest;
	algorithm::GraphicItemParameter m_region_toProcess;

	/// The intermediate image factory after color deconvolution
	std::shared_ptr<image::tile::Factory> m_colorDeconvolution_factory;

	/// The intermediate image factory after thresholding
	//std::shared_ptr<image::tile::Factory> m_threshold_factory;
	std::ofstream log_file;


};

} // namespace algorithm
} // namespace sedeen

#endif

