/*=========================================================================
 *
 *  Copyright (c) 2019 Sunnybrook Research Institute
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

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
#include <boost/filesystem.hpp>

using namespace boost::filesystem;


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
					///rgb_from_hex: Hematoxylin + Eosin
					HematoxylinPEosin,
					///rgb_from_hdx: Hematoxylin + DAB
					HematoxylinPDAB,
					///rgb_from_hed: Hematoxylin + Eosin + DAB
					HematoxylinPEosinPDAB,
					///Load from previously saved file
					LoadFromFile
				};

				/// Disply options for the output image
				enum DisplyOptions {
					STAIN1,
					STAIN2,
					STAIN3
				};

				/// Creates a colour deconvolution Kernel with selected 
				// color-deconvolution matrix
				//
				/// \param 
				/// 
				explicit ColorDeconvolution(Behavior behavior, DisplyOptions displyOption, double[9], double, const std::string& );

				virtual ~ColorDeconvolution();


				/// Set the display option of the kernel
				/// \param t
				/// The stain name
				/// \post
				///  is updated. If  is changed, update() is called to notify
				/// the observers
				void setDisplayOptions(DisplyOptions displyOption);

				/// Set the Stain combination matrices of the kernel
				/// \param t
				/// The Stain combination matrice name
				/// \post
				/// m_StainMatrice is updated. If m_StainMatrice is changed, update() is called to notify
				/// the observers.
				void SetStainMatrice(Behavior behavior, double[9], const std::string&);

			private:
				/// \cond INTERNAL

				virtual RawImage doProcessData(const RawImage &source);

				virtual const ColorSpace& doGetColorSpace() const;

				//RawImage separate_stains(const RawImage &source, const Eigen::Matrix<double,3,3>);
				/*Eigen::MatrixXd convertImageToMatrix(const RawImage &source);
				RawImage convertMatrixToImage(const Eigen::MatrixXd &source);*/

				RawImage separate_stains(const RawImage &source, double[9]);
				void computeMatrixInvers( double[9] );
				bool saveToCSVfile(const std::string& );
				void loadFromCSV( const std::string&, const std::string&);
				/*void getmeanRGBODfromROI(std::shared_ptr<tile::Factory> source,
					const Rect &region_of_interest,
					const Size &rescaled_resolution);*/
				///matrix conversion from H&E to RGB (original matrix from Ruifrok, 2001)
				// rows of matrix are stains, columns are color channels
				//const Eigen::Matrix<double, 3, 3, Eigen::DontAlign> m_rgb_from_HandE;
				ColorDeconvolution::Behavior m_Stain;
				ColorDeconvolution::DisplyOptions m_displyOption;	
				//std::vector<RawImage> m_outputImages;
				//std::vector<RawImage> m_binaryImages;
				double m_threshold;

				static const int NumOfStains =3;
				double log255;
				double MODx[3];
				double MODy[3];
				double MODz[3];
				double cosx[3];
				double cosy[3];
				double cosz[3];

				//std::ofstream log_file;
				int count;

        ColorSpace m_colorSpace;

				/// \endcond
			};

		} // namespace tile

		PATHCORE_IMAGE_API
			void getStainsComponents(std::shared_ptr<tile::Factory> source,
								const std::vector<std::shared_ptr<GraphicItemBase>> region_of_interests,
								const Size& rescaled_resolutions, double[9]);

		PATHCORE_IMAGE_API
			void getmeanRGBODfromROI(RawImage, double[3]);

		PATHCORE_IMAGE_API
		bool serachForfile( const path&, const std::string&, path& ); 

		PATHCORE_IMAGE_API
		void setImagePath( const std::string& ); 

	} // namespace image
} // namespace sedeen
#endif

