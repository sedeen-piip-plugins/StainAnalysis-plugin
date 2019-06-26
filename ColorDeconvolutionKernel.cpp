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

#include "ColorDeconvolutionKernel.h"
#include "StainVectorMath.h"

namespace sedeen {
	namespace image {
		namespace tile {
			ColorDeconvolution::ColorDeconvolution( DisplayOptions displayOption, 
                std::shared_ptr<StainProfile> theProfile, 
                bool applyThreshold, double threshold) :
                m_applyThreshold(applyThreshold),
				m_threshold(threshold),
				m_DisplayOption(displayOption),
                m_stainProfile(theProfile),
                m_colorSpace(ColorModel::RGBA, ChannelType::UInt8)
			{
			}//end constructor

            ColorDeconvolution::~ColorDeconvolution(void) {
            }//end destructor

			RawImage ColorDeconvolution::separateStains(const RawImage &source, double stainVec_matrix[9])
			{
                int scaleMax = 255;
                // initialize 3 output images
                sedeen::Size imageSize = source.size();
                std::vector<RawImage> colorImages;
                for (int i = 0; i < 3; i++) {
                    //Color images: adjust pixel value, assign
                    colorImages.push_back(RawImage(imageSize, ColorSpace(ColorModel::RGBA, ChannelType::UInt8)));
                    colorImages[i].fill(0);
                }
                //Get the number of stains in the profile
                //int numStains = m_stainProfile->GetNumberOfStainComponents();

                //The inverse can't be calculated if there is a row of zeros. Replace these values first.
                //If there is a row of zeros, replace it. Use the default replacement row.
                double noZeroRowsMatrix[9] = { 0.0 };
                StainVectorMath::ConvertZeroRowsToUnitary(stainVec_matrix, noZeroRowsMatrix);

                //Get the inverse of the noZeroRowsMatrix
                double inverse_matrix[9] = { 0.0 };
                StainVectorMath::Compute3x3MatrixInverse(noZeroRowsMatrix, inverse_matrix);

                //double normalized_inverse_matrix[9] = { 0.0 };
                //StainVectorMath::Make3x3MatrixUnitary(inverse_matrix, normalized_inverse_matrix);


//                //Open a file for temporary output
//#include <fstream>
//                std::fstream tempOut;
//                tempOut.open("D:\\mschumaker\\projects\\Sedeen\\testData\\output\\tempout.txt", std::fstream::out);
//
//                std::stringstream sv;
//                sv << "The stain vector matrix:\n";
//                sv << stainVec_matrix[0] << ", " << stainVec_matrix[1] << ", " << stainVec_matrix[2] << "\n";
//                sv << stainVec_matrix[3] << ", " << stainVec_matrix[4] << ", " << stainVec_matrix[5] << "\n";
//                sv << stainVec_matrix[6] << ", " << stainVec_matrix[7] << ", " << stainVec_matrix[8];
//                tempOut << sv.str() << std::endl;
//
//                std::stringstream sz;
//                sz << "The no-zero-rows matrix:\n";
//                sz << noZeroRowsMatrix[0] << ", " << noZeroRowsMatrix[1] << ", " << noZeroRowsMatrix[2] << "\n";
//                sz << noZeroRowsMatrix[3] << ", " << noZeroRowsMatrix[4] << ", " << noZeroRowsMatrix[5] << "\n";
//                sz << noZeroRowsMatrix[6] << ", " << noZeroRowsMatrix[7] << ", " << noZeroRowsMatrix[8];
//                tempOut << sz.str() << std::endl;
//
//                std::stringstream si;
//                si << "The inverse matrix:\n";
//                si << inverse_matrix[0] << ", " << inverse_matrix[1] << ", " << inverse_matrix[2] << "\n";
//                si << inverse_matrix[3] << ", " << inverse_matrix[4] << ", " << inverse_matrix[5] << "\n";
//                si << inverse_matrix[6] << ", " << inverse_matrix[7] << ", " << inverse_matrix[8];
//                tempOut << si.str() << std::endl;
//
//                std::stringstream sn;
//                sn << "The normalized inverse matrix:\n";
//                sn << normalized_inverse_matrix[0] << ", " << normalized_inverse_matrix[1] << ", " << normalized_inverse_matrix[2] << "\n";
//                sn << normalized_inverse_matrix[3] << ", " << normalized_inverse_matrix[4] << ", " << normalized_inverse_matrix[5] << "\n";
//                sn << normalized_inverse_matrix[6] << ", " << normalized_inverse_matrix[7] << ", " << normalized_inverse_matrix[8];
//                tempOut << sn.str() << std::endl;
//
//                //close the temporary output file
//                tempOut.close();


                //loop over all pixels in the given RawImage		
				int y = 0, x = 0;
                for (int j = 0; j < imageSize.width()*imageSize.height(); j++) {
					x = j%imageSize.width();
					y = j/imageSize.width();

					// log transform the RGB data
                    double R = source.at(x, y, 0).as<double>();
                    double G = source.at(x, y, 1).as<double>();
                    double B = source.at(x, y, 2).as<double>();
                    double pixelOD[3]; //index is color channel
                    pixelOD[0] = StainVectorMath::ConvertRGBtoOD(R);
                    pixelOD[1] = StainVectorMath::ConvertRGBtoOD(G);
                    pixelOD[2] = StainVectorMath::ConvertRGBtoOD(B);

                    //The resulting RGB values for the three images
                    double RGB_sep[9] = { 0.0 };
                    GetSeparateColorsForPixel(pixelOD, RGB_sep, stainVec_matrix, inverse_matrix);

                    for (int i = 0; i < 3; i++) { //i index is stain
                        colorImages.at(i).setValue(x, y, 0, static_cast<int>(RGB_sep[i * 3    ]));
                        colorImages.at(i).setValue(x, y, 1, static_cast<int>(RGB_sep[i * 3 + 1]));
                        colorImages.at(i).setValue(x, y, 2, static_cast<int>(RGB_sep[i * 3 + 2]));
                        colorImages.at(i).setValue(x, y, 3, scaleMax);
                    }
				}//end for each pixel

                //Return the requested stain image
				if( m_DisplayOption == DisplayOptions::STAIN1 ){
					return colorImages[0];
                }
				else if( m_DisplayOption == DisplayOptions::STAIN2 ){
					return colorImages[1];
                }
				else if( m_DisplayOption == DisplayOptions::STAIN3 ){
					return colorImages[2];
                }
                else {
                    return source;
                }
			}//end separateStains

            void ColorDeconvolution::GetSeparateColorsForPixel(double pixelOD[3], double RGB_sep[9], double stainVec_matrix[9], double inverse_matrix[9]) {
                //Determine how much of each stain is present at a pixel
                double stainSaturation[3] = { 0.0 }; //index is stain number
                StainVectorMath::Multiply3x3MatrixAndVector(inverse_matrix, pixelOD, stainSaturation);

                for (int i = 0; i < 3; i++) { //i index is stain
                    //Scale the stain's OD by the amount of stain at this pixel, get the RGB values
                    double OD_scaled[3]; //index is color channel

                    //Don't allow negative stain saturations
                    stainSaturation[i] = (stainSaturation[i] > 0.0) ? stainSaturation[i] : 0.0;

                    OD_scaled[0] = (stainSaturation[i]) * stainVec_matrix[i * 3];
                    OD_scaled[1] = (stainSaturation[i]) * stainVec_matrix[i * 3 + 1];
                    OD_scaled[2] = (stainSaturation[i]) * stainVec_matrix[i * 3 + 2];

                    double OD_sum = OD_scaled[0] + OD_scaled[1] + OD_scaled[2];
                    //Determine if the threshold should be applied to this stain's pixel value
                    bool isAboveThreshold;
                    if (m_applyThreshold) {
                        isAboveThreshold = (OD_sum > m_threshold) ? true : false;
                    }
                    else {
                        isAboveThreshold = true;
                    }

                    if (isAboveThreshold) {
                        RGB_sep[i * 3    ] = StainVectorMath::ConvertODtoRGB(OD_scaled[0]);
                        RGB_sep[i * 3 + 1] = StainVectorMath::ConvertODtoRGB(OD_scaled[1]);
                        RGB_sep[i * 3 + 2] = StainVectorMath::ConvertODtoRGB(OD_scaled[2]);
                    }
                    else {
                        RGB_sep[i * 3    ] = 0.0;
                        RGB_sep[i * 3 + 1] = 0.0;
                        RGB_sep[i * 3 + 2] = 0.0;
                    }
                }
            }//end GetSeparateColorsForPixel

            RawImage ColorDeconvolution::thresholdOnly(const RawImage &source) {
                int scaleMax = 255;
                // initialize 3 output images
                sedeen::Size imageSize = source.size();
                std::vector<RawImage> colorImages;
                for (int i = 0; i < 3; i++) {
                    //Color images: adjust pixel value, assign
                    colorImages.push_back(RawImage(imageSize, ColorSpace(ColorModel::RGBA, ChannelType::UInt8)));
                    colorImages[i].fill(0);
                }

                // translate ------------------				
                int y = 0, x = 0;
                for (int j = 0; j < imageSize.width()*imageSize.height(); j++) {
                    x = j % imageSize.width();
                    y = j / imageSize.width();

                    // log transform the RGB data
                    double R = source.at(x, y, 0).as<double>();
                    double G = source.at(x, y, 1).as<double>();
                    double B = source.at(x, y, 2).as<double>();
                    double pixelOD[3]; //index is color channel
                    pixelOD[0] = StainVectorMath::ConvertRGBtoOD(R);
                    pixelOD[1] = StainVectorMath::ConvertRGBtoOD(G);
                    pixelOD[2] = StainVectorMath::ConvertRGBtoOD(B);
                    //Get the total OD at the pixel
                    double OD_sum = pixelOD[0] + pixelOD[1] + pixelOD[2];

                    bool isAboveThreshold;
                    if (m_applyThreshold) {
                        isAboveThreshold = (OD_sum > m_threshold) ? true : false;
                    }
                    else {
                        isAboveThreshold = true;
                    }

                    //If the pixel was determined to be above threshold, add it to the images
                    for (int i = 0; i < 3; i++) {
                        if (isAboveThreshold) {
                            colorImages.at(i).setValue(x, y, 0, static_cast<int>(R));
                            colorImages.at(i).setValue(x, y, 1, static_cast<int>(G));
                            colorImages.at(i).setValue(x, y, 2, static_cast<int>(B));
                            colorImages.at(i).setValue(x, y, 3, scaleMax);
                        }
                        else {
                            colorImages.at(i).setValue(x, y, 0, 0);
                            colorImages.at(i).setValue(x, y, 1, 0);
                            colorImages.at(i).setValue(x, y, 2, 0);
                            colorImages.at(i).setValue(x, y, 3, scaleMax);
                        }
                    }
                }

                //Return the requested stain image
                if (m_DisplayOption == DisplayOptions::STAIN1) {
                    return colorImages[0];
                }
                else if (m_DisplayOption == DisplayOptions::STAIN2) {
                    return colorImages[1];
                }
                else if (m_DisplayOption == DisplayOptions::STAIN3) {
                    return colorImages[2];
                }
                else {
                    return source;
                }
            }//end thresholdOnly

			RawImage ColorDeconvolution::doProcessData(const RawImage &source)
			{
                double stainVec_matrix[9] = { 0.0 };
                //Fill stainVec_matrix with the stain vector profile values
                bool checkResult = m_stainProfile->GetNormalizedProfilesAsDoubleArray(stainVec_matrix);
                //Get the number of stains in the profile
                int numStains = m_stainProfile->GetNumberOfStainComponents();
                if (checkResult) {
                    //Is number of stains set to 1? Threshold only if so
                    if (numStains == 1) {
                        return thresholdOnly(source);
                    }
                    else if (numStains == 2 || numStains == 3) {
                        //Stain separation and thresholding
                        return separateStains(source, stainVec_matrix);
                    }
                    else {
                        return source;
                    }
                }
                else {
                    return source;
                }
			}//end doProcessData

			const ColorSpace& ColorDeconvolution::doGetColorSpace() const
			{
				return m_colorSpace;
			}

		} // namespace tile


		void getmeanRGBODfromROI(RawImage ROI, double rgbOD[3])
		{
			if( ROI.isNull() )
				return;
            //temporary array
            double tempOD[3] = { 0.0 };

			int imageSize = ROI.size().width()*ROI.size().height();
			double log255= log(255.0);
			int y = 0, x = 0;
			for (int i = 0; i < imageSize; i++){
				x = i%ROI.size().width();
                if (x == 0 && i != 0) {
                    y++;
                }
                //Convert RGB vals to optical density, sum over all pixels
                tempOD[0] = tempOD[0] + StainVectorMath::ConvertRGBtoOD(ROI.at(x, y, 0).as<double>());
                tempOD[1] = tempOD[1] + StainVectorMath::ConvertRGBtoOD(ROI.at(x, y, 1).as<double>());
                tempOD[2] = tempOD[2] + StainVectorMath::ConvertRGBtoOD(ROI.at(x, y, 2).as<double>());
			}
            //average of all pixels in region of interest
            rgbOD[0] = tempOD[0] / imageSize;
			rgbOD[1] = tempOD[1] / imageSize;
			rgbOD[2] = tempOD[2] / imageSize;
		}//end getmeanRGBODfromROI

		void getStainsComponents(std::shared_ptr<tile::Factory> source,
			const std::vector<std::shared_ptr<GraphicItemBase>> region_of_interests,
			const Size &rescaled_resolutions, double stainVec_matrix[9])
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
			for(size_t j=0; j < numberOfRegions; j++)
			{
				rgbOD[0]=0.0;
				rgbOD[1]=0.0;
				rgbOD[2]=0.0;
				Rect rect = containingRect(region_of_interests.at(j)->graphic());
				RawImage ROI = compositor.getImage(rect, Size(rect.width(), rect.height()));
				getmeanRGBODfromROI(ROI, rgbOD);

                stainVec_matrix[j * 3    ] = rgbOD[0];
                stainVec_matrix[j * 3 + 1] = rgbOD[1];
                stainVec_matrix[j * 3 + 2] = rgbOD[2];
			}
		}//end getStainsComponents

	} // namespace image
} // namespace sedeen