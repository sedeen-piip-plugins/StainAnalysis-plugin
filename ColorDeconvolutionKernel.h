/*=============================================================================
 *
 *  Copyright (c) 2021 Sunnybrook Research Institute
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

#ifndef SEDEEN_SRC_IMAGE_FILTER_KERNELS_COLORDECONVOLUTION_H
#define SEDEEN_SRC_IMAGE_FILTER_KERNELS_COLORDECONVOLUTION_H

#include "Global.h"
#include "Geometry.h"
#include "Image.h"
#include "image/filter/Kernel.h"
#include "global/ColorSpace.h"

#include <cstdio>
#include <sstream>
#include <memory>
#include <filesystem> //Requires C++17

//Plugin includes
#include "StainProfile.h"

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
        explicit ColorDeconvolution(DisplayOptions displayOption, std::shared_ptr<StainProfile>, 
            bool applyThreshold, double threshold, bool stainQuantityOnly = false);

		virtual ~ColorDeconvolution();

	private:
		/// \cond INTERNAL

		virtual RawImage doProcessData(const RawImage &source);
		virtual const ColorSpace& doGetColorSpace() const;
        ///Given an input RGBA image and stain vectors, output either an RGBA or a Grayscale image
		RawImage separateStains(const RawImage &source, double (&stainVec)[9]);

        ///Apply threshold to source image, output RGB or averaged grayscale if stainQuantityOnly is false
        RawImage thresholdOnly(const RawImage &source);


        ///Arguments are: the 3 OD values, the 9-element RGB output, the 3-element stain quantity output, the stain vector matrix, and the inverse of the matrix
        void GetSeparateColorsForPixel(double(&pixelOD)[3], double(&RGB_sep)[9],
            double(&quant)[3], double(&stainVec_matrix)[9], double(&inverse_matrix)[9]);

        ///Original argument list (overload)
        void GetSeparateColorsForPixel(double (&pixelOD)[3], double (&RGB_sep)[9], 
            double (&stainVec_matrix)[9], double (&inverse_matrix)[9]);

        ///Get the boolean value of m_grayscaleQuanityOnly
        const bool& GetGrayscaleQuantityOnly() const { return m_grayscaleQuantityOnly; }
        ///Set the output color space
        void SetOutputColorSpace(const ColorSpace& os) { m_outputColorSpace = os; }

		// rows of matrix are stains, columns are color channels
		ColorDeconvolution::DisplayOptions m_DisplayOption;	
        bool m_grayscaleQuantityOnly;
        bool m_applyThreshold;
		double m_threshold;

        //Normalization for grayscale stain quantities: -log_10(1/255) is 2.40, so set 2.55 to be channel 255 = norm factor 100
        const double m_grayscaleNormFactor;

        ///ColorSpace of the source (input) image
        ///ColorSpace of the output image
        ColorSpace m_outputColorSpace;
        std::shared_ptr<StainProfile> m_stainProfile;
		/// \endcond
	};

} // namespace tile
} // namespace image
} // namespace sedeen
#endif