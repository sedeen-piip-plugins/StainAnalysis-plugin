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
			ColorDeconvolution::ColorDeconvolution( Behavior behavior, DisplayOptions displayOption, 
                std::shared_ptr<StainProfile> theProfile, 
                bool applyThreshold, double threshold) :
                m_applyThreshold(applyThreshold),
				m_threshold(threshold),
				m_behaviorType(behavior),
				m_DisplayOption(displayOption),
                m_stainProfile(theProfile),
                m_colorSpace(ColorModel::RGBA, ChannelType::UInt8),
                m_fullColorImages(true) //return full color images (true) or binary images (false)
			{
			}//end constructor

            ColorDeconvolution::~ColorDeconvolution(void) {
            }//end destructor

			void ColorDeconvolution::computeMatrixInverse( double inputMat[9], double inversionMat[9] )
			{
				double leng, A, V, C;
				double  len[3];
                double MODx[3];
                double MODy[3];
                double MODz[3];
                double cosx[3];
                double cosy[3];
                double cosz[3];

                for (int i = 0; i < 3; i++) {
                    MODx[i] = inputMat[i * 3    ];
                    MODy[i] = inputMat[i * 3 + 1];
                    MODz[i] = inputMat[i * 3 + 2];
                }

                for (int i = 0; i < 3; i++) {
                    //normalize vector length
                    cosx[i] = cosy[i] = cosz[i] = 0.0;
                    len[i] = std::sqrt(MODx[i] * MODx[i] + MODy[i] * MODy[i] + MODz[i] * MODz[i]);
                    if (len[i] != 0.0) {
                        cosx[i] = MODx[i] / len[i];
                        cosy[i] = MODy[i] / len[i];
                        cosz[i] = MODz[i] / len[i];
                    }
                }

				// translation matrix
                if (cosx[1] == 0.0) { //2nd colour is unspecified
                    if (cosy[1] == 0.0) {
                        if (cosz[1] == 0.0) {
                            cosx[1] = cosz[0];
                            cosy[1] = cosx[0];
                            cosz[1] = cosy[0];
                        }
                    }
                }

				if (cosx[2] == 0.0) { // 3rd colour is unspecified
					if (cosy[2] == 0.0) {
						if (cosz[2] == 0.0) {
                            cosx[2]= cosy[0]*cosz[1] - cosy[1]*cosz[0];
                            cosy[2]= cosz[0]*cosx[1] - cosz[1]*cosx[0];
                            cosz[2]= cosx[0]*cosy[1] - cosx[1]*cosy[0];
						}
					}
				}

				leng=std::sqrt(cosx[2]*cosx[2] + cosy[2]*cosy[2] + cosz[2]*cosz[2]);

				cosx[2]= cosx[2]/leng;
				cosy[2]= cosy[2]/leng;
				cosz[2]= cosz[2]/leng;

				for (int i=0; i<3; i++) {
					if (cosx[i] == 0.0) cosx[i] = 0.001;
					if (cosy[i] == 0.0) cosy[i] = 0.001;
					if (cosz[i] == 0.0) cosz[i] = 0.001;
				}

				//matrix inversion
				A = cosy[1] - cosx[1] * cosy[0] / cosx[0];
				//if(A==0) A=0.001;
				V = cosz[1] - cosx[1] * cosz[0] / cosx[0];
				C = cosz[2] - cosy[2] * V/A + cosx[2] * (V/A * cosy[0] / cosx[0] - cosz[0] / cosx[0]);
				//if(C==0) C=0.001;
				inversionMat[2] = (-cosx[2] / cosx[0] - cosx[2] / A * cosx[1] / cosx[0] * cosy[0] / cosx[0] + cosy[2] / A * cosx[1] / cosx[0]) / C;
				inversionMat[1] = -inversionMat[2] * V / A - cosx[1] / (cosx[0] * A);
				inversionMat[0] = 1.0 / cosx[0] - inversionMat[1] * cosy[0] / cosx[0] - inversionMat[2] * cosz[0] / cosx[0];
				inversionMat[5] = (-cosy[2] / A + cosx[2] / A * cosy[0] / cosx[0]) / C;
				inversionMat[4] = -inversionMat[5] * V / A + 1.0 / A;
				inversionMat[3] = -inversionMat[4] * cosy[0] / cosx[0] - inversionMat[5] * cosz[0] / cosx[0];
				inversionMat[8] = 1.0 / C;
				inversionMat[7] = -inversionMat[8] * V / A;
				inversionMat[6] = -inversionMat[7] * cosy[0] / cosx[0] - inversionMat[8] * cosz[0] / cosx[0];
			}//end computeMatrixInverse

			RawImage ColorDeconvolution::separate_stains(const RawImage &source, double stainVec_matrix[9])
			{
                int scaleMax = 255;
                double inverse_matrix[9] = { 0.0 };
                //Get the inverse of the matrix
                computeMatrixInverse(stainVec_matrix, inverse_matrix);

				// initialize 3 output colour images
				sedeen::Size imageSize = source.size();
				std::vector<RawImage> binaryImages;
                std::vector<RawImage> colorImages;
				for (int i=0; i<3; i++){
                    //The binary images: all or nothing
					binaryImages.push_back(RawImage(imageSize, ColorSpace(ColorModel::RGBA, ChannelType::UInt8)));
					binaryImages[i].fill(0);
                    //Color images: adjust pixel value, assign
                    colorImages.push_back(RawImage(imageSize, ColorSpace(ColorModel::RGBA, ChannelType::UInt8)));
                    colorImages[i].fill(0);
                }

                //Get colors to assign to binary images from stainVec_matrix
                std::vector<std::array<double, 3>> binaryStainColors;
                for (int i = 0; i < 3; i++) {
                    std::array<double, 3> rgb;
                    rgb[0] = StainVectorMath::convertODtoRGB(stainVec_matrix[i * 3    ]);
                    rgb[1] = StainVectorMath::convertODtoRGB(stainVec_matrix[i * 3 + 1]);
                    rgb[2] = StainVectorMath::convertODtoRGB(stainVec_matrix[i * 3 + 2]);
                    //brighten the color (MaximizeArray returns values from 0 to 1
                    std::array<double, 3> rgb_max = StainVectorMath::MaximizeArray(rgb);
                    for (auto p = rgb_max.begin(); p != rgb_max.end(); ++p) { *p = *p * static_cast<double>(scaleMax); }
                    //push to the vector
                    binaryStainColors.push_back(rgb);
                }

				// translate ------------------				
				int y = 0, x = 0;
				for (int j=0;j<imageSize.width()*imageSize.height();j++){
					x = j%imageSize.width();
					y = j/imageSize.width();

					// log transform the RGB data
                    double R = source.at(x, y, 0).as<double>();
                    double G = source.at(x, y, 1).as<double>();
                    double B = source.at(x, y, 2).as<double>();

					double OD_R = StainVectorMath::convertRGBtoOD(R);
					double OD_G = StainVectorMath::convertRGBtoOD(G);
					double OD_B = StainVectorMath::convertRGBtoOD(B);

					for (int i=0; i<3; i++){
						// rescale to match original paper values
                        double OD_Rscaled = OD_R * inverse_matrix[i * 3    ];
                        double OD_Gscaled = OD_G * inverse_matrix[i * 3 + 1];
                        double OD_Bscaled = OD_B * inverse_matrix[i * 3 + 2];

                        //Create display colors
                        std::array<double, 3> RGB_sep;
                        RGB_sep[0] = StainVectorMath::convertODtoRGB(OD_Rscaled);
                        RGB_sep[1] = StainVectorMath::convertODtoRGB(OD_Gscaled);
                        RGB_sep[2] = StainVectorMath::convertODtoRGB(OD_Bscaled);

                        unsigned char binaryValue;
                        if (m_applyThreshold) {
                            binaryValue = ((OD_Rscaled + OD_Gscaled + OD_Bscaled) > m_threshold) ? 1 : 0;
                        }
                        else {
                            binaryValue = 1;
                        }

                        //If the pixel was determined to be above threshold, add it to the image
                        if (binaryValue) {
                            binaryImages.at(i).setValue(x, y, 0, static_cast<int>((binaryStainColors[i])[0]));
                            binaryImages.at(i).setValue(x, y, 1, static_cast<int>((binaryStainColors[i])[1]));
                            binaryImages.at(i).setValue(x, y, 2, static_cast<int>((binaryStainColors[i])[2]));
                            binaryImages.at(i).setValue(x, y, 3, scaleMax);
                            colorImages.at(i).setValue(x, y, 0, static_cast<int>(RGB_sep[0]));
                            colorImages.at(i).setValue(x, y, 1, static_cast<int>(RGB_sep[1]));
                            colorImages.at(i).setValue(x, y, 2, static_cast<int>(RGB_sep[2]));
                            colorImages.at(i).setValue(x, y, 3, scaleMax);
                        }
                        else {
                            binaryImages.at(i).setValue(x, y, 0, 0);
                            binaryImages.at(i).setValue(x, y, 1, 0);
                            binaryImages.at(i).setValue(x, y, 2, 0);
                            binaryImages.at(i).setValue(x, y, 3, scaleMax);
                            colorImages.at(i).setValue(x, y, 0, 0);
                            colorImages.at(i).setValue(x, y, 1, 0);
                            colorImages.at(i).setValue(x, y, 2, 0);
                            colorImages.at(i).setValue(x, y, 3, scaleMax);
                        }
					}
				}

                //Get the value of the member variable m_fullColorImages to decide return type
                bool returnFullColor = m_fullColorImages;
				if( m_DisplayOption == DisplayOptions::STAIN1 ){
					return returnFullColor ? colorImages[0] : binaryImages[0];
                }
				else if( m_DisplayOption == DisplayOptions::STAIN2 ){
					return returnFullColor ? colorImages[1] : binaryImages[1];
                }
				else if( m_DisplayOption == DisplayOptions::STAIN3 ){
					return returnFullColor ? colorImages[2] : binaryImages[2];
                }
                else {
                    return source;
                }
			}//end separate_stains

			RawImage ColorDeconvolution::doProcessData(const RawImage &source)
			{
                double stainVec_matrix[9] = { 0.0 };
                //Fill stainVec_matrix with the stain vector profile values
                bool checkResult = m_stainProfile->GetNormalizedProfilesAsDoubleArray(stainVec_matrix);
                if (checkResult) {
                    //Stain separation and thresholding
                    return separate_stains(source, stainVec_matrix);
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
                tempOD[0] = tempOD[0] + StainVectorMath::convertRGBtoOD(ROI.at(x, y, 0).as<double>());
                tempOD[1] = tempOD[1] + StainVectorMath::convertRGBtoOD(ROI.at(x, y, 1).as<double>());
                tempOD[2] = tempOD[2] + StainVectorMath::convertRGBtoOD(ROI.at(x, y, 2).as<double>());
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