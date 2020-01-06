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

#include "ColorDeconvolutionKernel.h"
#include "StainRGBFromROI.h"
#include "StainVectorMath.h"
#include "ODConversion.h"

namespace sedeen {
    namespace image {


        void getmeanRGBODfromROI(RawImage ROI, double (&rgbOD)[3])
        {
            if (ROI.isNull())
                return;
            //temporary array
            double tempOD[3] = { 0.0 };

            //Perform fast color to OD conversion using a lookup table
            std::shared_ptr<ODConversion> converter = std::make_shared<ODConversion>();

            int imageSize = ROI.size().width()*ROI.size().height();
            double log255 = log(255.0);
            int y = 0, x = 0;
            for (int i = 0; i < imageSize; i++) {
                x = i % ROI.size().width();
                if (x == 0 && i != 0) {
                    y++;
                }
                //Convert RGB vals to optical density, sum over all pixels
                tempOD[0] = tempOD[0] + converter->LookupRGBtoOD(ROI.at(x, y, 0).as<int>());
                tempOD[1] = tempOD[1] + converter->LookupRGBtoOD(ROI.at(x, y, 1).as<int>());
                tempOD[2] = tempOD[2] + converter->LookupRGBtoOD(ROI.at(x, y, 2).as<int>());
            }
            //average of all pixels in region of interest
            rgbOD[0] = tempOD[0] / imageSize;
            rgbOD[1] = tempOD[1] / imageSize;
            rgbOD[2] = tempOD[2] / imageSize;
        }//end getmeanRGBODfromROI

        void getStainsComponents(std::shared_ptr<tile::Factory> source,
            const std::vector<std::shared_ptr<GraphicItemBase>> region_of_interests,
            const Size &rescaled_resolutions, double (&stainVec_matrix)[9])
        {
            for (int i = 0; i < 9; i++) {
                stainVec_matrix[i] = 0;
            }

            //Get the length of the region_of_interests vector, though we want at most 3
            size_t numberOfRegions = region_of_interests.size();
            numberOfRegions = (numberOfRegions > 3) ? 3 : numberOfRegions;

            // Get image from the output factory
            auto compositor = image::tile::Compositor(source);
            double  rgbOD[3];
            for (size_t j = 0; j < numberOfRegions; j++)
            {
                rgbOD[0] = 0.0;
                rgbOD[1] = 0.0;
                rgbOD[2] = 0.0;
                Rect rect = containingRect(region_of_interests.at(j)->graphic());
                RawImage ROI = compositor.getImage(rect, Size(rect.width(), rect.height()));
                getmeanRGBODfromROI(ROI, rgbOD);

                stainVec_matrix[j * 3] = rgbOD[0];
                stainVec_matrix[j * 3 + 1] = rgbOD[1];
                stainVec_matrix[j * 3 + 2] = rgbOD[2];
            }
        }//end getStainsComponents


    } // namespace image
} // namespace sedeen
