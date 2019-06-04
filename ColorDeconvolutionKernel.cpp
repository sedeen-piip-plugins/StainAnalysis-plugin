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

namespace sedeen {
	namespace image {
		namespace tile {
			ColorDeconvolution::ColorDeconvolution( Behavior behavior, DisplayOptions displayOption, 
                std::shared_ptr<StainProfile> theProfile, 
                bool applyThreshold, double threshold) :
                m_applyThreshold(applyThreshold),
				m_threshold(threshold),
				m_behaviorType(behavior),
				//log255(std::log(255.0)),
				m_DisplayOption(displayOption),
                m_stainProfile(theProfile),
                m_colorSpace(ColorModel::RGBA, ChannelType::UInt8)
			{
                //Get the stain profile values as a 9-element array
                //double conv_matrix[9];
                //theProfile->GetProfileAsDoubleArray(conv_matrix);
                //SetStainMatrix(m_behaviorType, conv_matrix);
				//setDisplayOptions(displayOption);
			}//end constructor


            ColorDeconvolution::~ColorDeconvolution(void) {
            }//end destructor

			//void ColorDeconvolution::setDisplayOptions(DisplayOptions displayOption)
			//{
   //             if (displayOption == DisplayOptions::STAIN1) {
   //                 m_DisplayOption = DisplayOptions::STAIN1;
   //             }

   //             if (displayOption == DisplayOptions::STAIN2) {
   //                 m_DisplayOption = DisplayOptions::STAIN2;
   //             }

   //             if (displayOption == DisplayOptions::STAIN3) {
   //                 m_DisplayOption = DisplayOptions::STAIN3;
   //             }
			//}//end setDisplayOptions

			//void ColorDeconvolution::SetStainMatrix(Behavior behavior, double conv_matrix[9])
			//{
   //             //std::filesystem::path dir_path = path_to_root +"sedeen/";
			//	if(behavior == Behavior::RegionOfInterest)
			//	{
			//		for (int i=0; i < 3; i++)
			//		{
			//			m_MODx[i]= conv_matrix[i*3];
			//			m_MODy[i]= conv_matrix[i*3+1];
			//			m_MODz[i]= conv_matrix[i*3+2];
			//		}

   //                 m_stainProfile->SetProfileFromDoubleArray(conv_matrix);
			//	}

			//	if(behavior == Behavior::LoadFromProfile){
   //                 //Fill conv_matrix with the stain profile values, then assign to member MOD arrays
   //                 bool checkResult = m_stainProfile->GetProfileAsDoubleArray(conv_matrix);
   //                 //If checkResult is false, assign zeros to everything
   //                 m_MODx[0] = checkResult ? conv_matrix[0] : 0.0;
			//		m_MODy[0] = checkResult ? conv_matrix[1] : 0.0;
			//		m_MODz[0] = checkResult ? conv_matrix[2] : 0.0;

   //                 m_MODx[1] = checkResult ? conv_matrix[3] : 0.0;
			//		m_MODy[1] = checkResult ? conv_matrix[4] : 0.0;
			//		m_MODz[1] = checkResult ? conv_matrix[5] : 0.0;

   //                 m_MODx[2] = checkResult ? conv_matrix[6] : 0.0;
			//		m_MODy[2] = checkResult ? conv_matrix[7] : 0.0;
			//		m_MODz[2] = checkResult ? conv_matrix[8] : 0.0;
			//	}
			//}//end SetStainMatrix

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


			RawImage ColorDeconvolution::separate_stains(const RawImage &source, double conv_matrix[9])
			{
				// initialize 3 output colour images
				sedeen::Size imageSize = source.size();
				std::vector<RawImage> binaryImages;
                //std::vector<RawImage> fullColorImages;
				for (int i=0; i<3; i++){
                    //The binary images: all or nothing
					binaryImages.push_back( RawImage(imageSize, ColorSpace(ColorModel::RGBA, ChannelType::UInt8)) );
					binaryImages[i].fill(0);
				}
				
				// translate ------------------				
				int y=0, x =0;
				for (int j=0;j<imageSize.width()*imageSize.height();j++){

					x = j%imageSize.width();
					y = j/imageSize.width();
					/*if( x == 0 && j!=0)
						y++;*/

					// log transform the RGB data
                    double R = source.at(x, y, 0).as<double>();
                    double G = source.at(x, y, 1).as<double>();
                    double B = source.at(x, y, 2).as<double>();
					double Rlog = -((255.0*log((R+1)/255.0))/log(255));
					double Glog = -((255.0*log((G+1)/255.0))/log(255));
					double Blog = -((255.0*log((B+1)/255.0))/log(255));

					double Rlogn = -log( (R/255.0 + 0.001)/1.001 );
					double Glogn = -log( (G/255.0 + 0.001)/1.001 );
					double Blogn = -log( (B/255.0 + 0.001)/1.001 );

					for (int i=0; i<3; i++){
						// rescale to match original paper values
						double Rscaled = Rlog * conv_matrix[i*3];
						double Gscaled = Glog * conv_matrix[i*3+1];
						double Bscaled = Blog * conv_matrix[i*3+2];

						double Rlogscaled = Rlogn * conv_matrix[i*3];
						double Glogscaled = Glogn * conv_matrix[i*3+1];
						double Blogscaled = Blogn * conv_matrix[i*3+2];

                        //Create display colors
                        int Rvector = static_cast<int>(255.0 * conv_matrix[i * 3]);
                        int Gvector = static_cast<int>(255.0 * conv_matrix[i * 3 + 1]);
                        int Bvector = static_cast<int>(255.0 * conv_matrix[i * 3 + 2]);

                        //Need a condition around this using m_applyThreshold
                        unsigned char binaryValue;
                        if (m_applyThreshold) {
                            binaryValue = ((Rlogscaled + Glogscaled + Blogscaled) > m_threshold) ? 255 : 0; 
                        }
                        else {
                            binaryValue = ((Rlogscaled + Glogscaled + Blogscaled) > 0.0001) ? 255 : 0;
                        }

						//binaryImages.at(i).setValue(x, y, 0, binaryValue);
						//binaryImages.at(i).setValue(x, y, 1, 0);
						//binaryImages.at(i).setValue(x, y, 2, 0);
						//binaryImages.at(i).setValue(x, y, 3, 255);
                        //If the pixel was determined to be above threshold, add it to the image
                        if (binaryValue) {
                            binaryImages.at(i).setValue(x, y, 0, Rvector);
                            binaryImages.at(i).setValue(x, y, 1, Gvector);
                            binaryImages.at(i).setValue(x, y, 2, Bvector);
                            binaryImages.at(i).setValue(x, y, 3, 255);
                        }
                        else {
                            binaryImages.at(i).setValue(x, y, 0, 0);
                            binaryImages.at(i).setValue(x, y, 1, 0);
                            binaryImages.at(i).setValue(x, y, 2, 0);
                            binaryImages.at(i).setValue(x, y, 3, 255);
                        }
					}
				}

				if( m_DisplayOption == DisplayOptions::STAIN1 ){					
					return binaryImages[0];
                }
				else if( m_DisplayOption == DisplayOptions::STAIN2 ){
					return binaryImages[1];
                }
				else if( m_DisplayOption == DisplayOptions::STAIN3 ){
					return binaryImages[2];
                }
                else {
                    return source;
                }
			}//end separate_stains

			RawImage ColorDeconvolution::doProcessData(const RawImage &source)
			{
                double stainVec_matrix[9] = { 0.0 };
                double inverse_matrix[9] = { 0.0 };

                //Fill stainVec_matrix with the stain vector profile values
                bool checkResult = m_stainProfile->GetNormalizedProfilesAsDoubleArray(stainVec_matrix);
                //Get the inverse of the matrix
                computeMatrixInverse(stainVec_matrix, inverse_matrix );

				//Stain Separation and combination
				return separate_stains( source, inverse_matrix );
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

			int imageSize = ROI.size().width()*ROI.size().height();
			double log255= log(255.0);
			int y = 0, x = 0;
			for (int i = 0; i < imageSize ; i++){
				x = i%ROI.size().width();
                if (x == 0 && i != 0) {
                    y++;
                }

				// rescale to match original paper values
				rgbOD[0] = rgbOD[0] + (-((255.0*log(((ROI.at(x, y, 0).as<double>() +1)/255.0))/log255)));
				rgbOD[1] = rgbOD[1] + (-((255.0*log(((ROI.at(x, y, 1).as<double>() +1)/255.0))/log255)));
				rgbOD[2] = rgbOD[2] + (-((255.0*log(((ROI.at(x, y, 2).as<double>() +1)/255.0))/log255)));
			}
			rgbOD[0] = rgbOD[0] / imageSize;
			rgbOD[1] = rgbOD[1] / imageSize;
			rgbOD[2] = rgbOD[2] / imageSize;
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