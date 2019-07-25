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

#ifndef DPTK_SRC_IMAGE_FILTER_KERNELS_COLORDECONVOLUTION_H
#define DPTK_SRC_IMAGE_FILTER_KERNELS_COLORDECONVOLUTION_H


#include <stdio.h>

//#include "image/filter/Kernel.h"
//#include "image\iterator\Iterator.h"
//#include "image\tile\Compositor.h"
#include "Global.h"
#include "Geometry.h"
#include "Image.h"

#include <fstream>
#include <sstream>
#include <filesystem> //Requires C++17
#include <memory>

//Plugin includes
#include "StainProfile.h"

//#include "Eigen\Dense"

namespace sedeen {

	namespace image {
		namespace tile {	

			/// \ingroup algorithm_kernels
			/// Color Deconvolution
			/// Apply colour deconvolution method described in:
			//Ruifrok AC, Johnston DA. Quantification of histochemical
			//staining by color deconvolution. Analytical & Quantitative
			//Cytology & Histology 2001; 23: 291-299
			class PATHCORE_IMAGE_API ColorDeconvolution : public Kernel {
			public:
				/// Stain combination matrices
				enum Behavior {
					///User selected ROI
					RegionOfInterest,
                    ///Load from the StainProfile
					LoadFromProfile
				};

				/// Display options for the output image
				enum DisplayOptions {
					STAIN1,
					STAIN2,
					STAIN3
				};

				/// Creates a colour deconvolution Kernel with selected 
				// color-deconvolution matrix
				//
				/// \param 
				/// 
                explicit ColorDeconvolution(Behavior behavior, DisplayOptions displayOption, std::shared_ptr<StainProfile>, bool, double);

				virtual ~ColorDeconvolution();


				/// Set the display option of the kernel
				/// \param t
				/// The stain name
				/// \post
				///  is updated. If  is changed, update() is called to notify
				/// the observers
				//void setDisplayOptions(DisplayOptions displayOption);

				/// Set the Stain combination matrices of the kernel
				/// \param t
				/// The Stain combination Matrix name
				/// \post
				/// m_StainMatrix is updated. If m_StainMatrix is changed, update() is called to notify
				/// the observers.
                //void SetStainMatrix(Behavior behavior, double[9]); // , const std::string&);

			private:
				/// \cond INTERNAL

				virtual RawImage doProcessData(const RawImage &source);

				virtual const ColorSpace& doGetColorSpace() const;

				//RawImage separate_stains(const RawImage &source, const Eigen::Matrix<double,3,3>);
				/*Eigen::MatrixXd convertImageToMatrix(const RawImage &source);
				RawImage convertMatrixToImage(const Eigen::MatrixXd &source);*/

				RawImage separate_stains(const RawImage &source, double[9]);
				void computeMatrixInverse( double[9], double[9] );
				///matrix conversion from H&E to RGB (original matrix from Ruifrok, 2001)
				// rows of matrix are stains, columns are color channels
				//const Eigen::Matrix<double, 3, 3, Eigen::DontAlign> m_rgb_from_HandE;
				ColorDeconvolution::Behavior m_behaviorType;
				ColorDeconvolution::DisplayOptions m_DisplayOption;	
				//std::vector<RawImage> m_outputImages;
				//std::vector<RawImage> m_binaryImages;
                bool m_applyThreshold;
				double m_threshold;

				//static const int NumOfStains = 3;
				//double log255;
				//double m_MODx[3];
				//double m_MODy[3];
				//double m_MODz[3];
				//double m_cosx[3];
				//double m_cosy[3];
				//double m_cosz[3];

				//std::ofstream log_file;
				//int count;

                ColorSpace m_colorSpace;

                std::shared_ptr<StainProfile> m_stainProfile;

				/// \endcond
			};

		} // namespace tile

		PATHCORE_IMAGE_API
			void getStainsComponents(std::shared_ptr<tile::Factory> source,
								const std::vector<std::shared_ptr<GraphicItemBase>> regions_of_interest,
								const Size& rescaled_resolutions, double[9]);

		PATHCORE_IMAGE_API
			void getmeanRGBODfromROI(RawImage, double[3]);

	} // namespace image
} // namespace sedeen
#endif

