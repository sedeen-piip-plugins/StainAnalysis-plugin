/*=========================================================================
 *
 *  Copyright (c) 2019 Sunnybrook Research Institute
 *
 *  License terms pending.
 *
 *=========================================================================*/

#include "ColorDeconvolutionKernel.h"

namespace sedeen {
	namespace image {
		namespace tile {
			ColorDeconvolution::ColorDeconvolution( Behavior behavior, DisplayOptions displayOption, std::shared_ptr<StainProfile> theProfile, 
                bool applyThreshold, double threshold) : /// const std::string& path_to_root) :
                m_applyThreshold(applyThreshold),
				m_threshold(threshold),
				count(0),
				m_behaviorType(behavior),
				log255(std::log(255.0)),
				m_DisplayOption(displayOption),
                m_stainProfile(theProfile),
                m_colorSpace(ColorModel::RGBA, ChannelType::UInt8)
			{
                //Get the stain profile values as a 9-element array
                double conv_matrix[9];
                theProfile->GetProfileAsDoubleArray(conv_matrix);
                SetStainMatrix(m_behaviorType, conv_matrix); // , path_to_root);
				setDisplayOptions(displayOption);
			}//end constructor


            ColorDeconvolution::~ColorDeconvolution(void) {
            }//end destructor

			void ColorDeconvolution::setDisplayOptions(DisplayOptions displayOption)
			{
                if (displayOption == DisplayOptions::STAIN1) {
                    m_DisplayOption = DisplayOptions::STAIN1;
                }

                if (displayOption == DisplayOptions::STAIN2) {
                    m_DisplayOption = DisplayOptions::STAIN2;
                }

                if (displayOption == DisplayOptions::STAIN3) {
                    m_DisplayOption = DisplayOptions::STAIN3;
                }
			}//end setDisplayOptions

			void ColorDeconvolution::SetStainMatrix(Behavior behavior, double conv_matrix[9]) //, const std::string& path_to_root)
			{
                //std::filesystem::path dir_path = path_to_root +"sedeen/";
				if(behavior == Behavior::RegionOfInterest)
				{
					for (int i=0; i < 3; i++)
					{
						m_MODx[i]= conv_matrix[i*3];
						m_MODy[i]= conv_matrix[i*3+1];
						m_MODz[i]= conv_matrix[i*3+2];
					}

                    //std::filesystem::path dir_path = path_to_root +"sedeen/";
					//std::string file_name = "StainsFile.csv";

					/*path dir_path = "C:/sedeen/";
					std::string file_name = "StainsFile.csv";*/
					//try
					//{
					//	//&& boost::filesystem::is_directory(dir_path)
					//	if( !std::filesystem::exists( dir_path ) )
					//	{
					//		if(std::filesystem::create_directory( dir_path ))
					//		{
					//			std::string fn = dir_path.string() + file_name;
					//			saveToCSVfile(fn);
					//		}
                    //
					//	}
					//	else
					//	{
					//		std::string fn = dir_path.string() + file_name;
					//		saveToCSVfile(fn);
					//	}
                    //
					//}
					//catch(std::filesystem::filesystem_error const & e)
					//{
                    //
					//}

                    m_stainProfile->SetProfileFromDoubleArray(conv_matrix);
				}

				if(behavior == Behavior::LoadFromProfile){
					//std::filesystem::path dir_path = "C:/sedeen/";
                    //std::filesystem::path dir_path = path_to_root +"sedeen/";
					//std::string file_name = "DefaultStainsFile.csv";
                    //std::filesystem::path path_found = std::filesystem::path();

					//bool thereIsAny = image::serachForfile( dir_path, file_name, path_found );
					//if(thereIsAny)
					//{
					//	std::string fn = dir_path.string() + file_name;
					//	loadFromCSV(fn, "HematoxylinPEosin");
					//}
					//else
					//{

                    //Fill conv_matrix with the stain profile values, then assign to member MOD arrays
                    bool checkResult = m_stainProfile->GetProfileAsDoubleArray(conv_matrix);
                    //If checkResult is false, assign zeros to everything
                    m_MODx[0] = checkResult ? conv_matrix[0] : 0.0;
					m_MODy[0] = checkResult ? conv_matrix[1] : 0.0;
					m_MODz[0] = checkResult ? conv_matrix[2] : 0.0;

                    m_MODx[1] = checkResult ? conv_matrix[3] : 0.0;
					m_MODy[1] = checkResult ? conv_matrix[4] : 0.0;
					m_MODz[1] = checkResult ? conv_matrix[5] : 0.0;

                    m_MODx[2] = checkResult ? conv_matrix[6] : 0.0;
					m_MODy[2] = checkResult ? conv_matrix[7] : 0.0;
					m_MODz[2] = checkResult ? conv_matrix[8] : 0.0;
					//}
				}
			}//end SetStainMatrix

			void ColorDeconvolution::computeMatrixInverse( double inversionMat[9] )
			{
				double leng, A, V, C;
				double  len[NumOfStains];
				//double  q[9];

				for (int i=0; i<NumOfStains; i++){
					//normalise vector length
					m_cosx[i]=m_cosy[i]=m_cosz[i]=0.0;
					len[i]=std::sqrt(m_MODx[i]*m_MODx[i] + m_MODy[i]*m_MODy[i] + m_MODz[i]*m_MODz[i]);
					if (len[i] != 0.0){
						m_cosx[i]= m_MODx[i]/len[i];
						m_cosy[i]= m_MODy[i]/len[i];
						m_cosz[i]= m_MODz[i]/len[i];
					}
				}

				// translation matrix
				if (m_cosx[1]==0.0){ //2nd colour is unspecified
					if (m_cosy[1]==0.0){
						if (m_cosz[1]==0.0){
							m_cosx[1]=m_cosz[0];
							m_cosy[1]=m_cosx[0];
							m_cosz[1]=m_cosy[0];
						}
					}
				}


				if (m_cosx[2]==0.0){ // 3rd colour is unspecified
					if (m_cosy[2]==0.0){
						if (m_cosz[2]==0.0){

							m_cosx[2]= m_cosy[0]*m_cosz[1] - m_cosy[1]*m_cosz[0];
							m_cosy[2]= m_cosz[0]*m_cosx[1] - m_cosz[1]*m_cosx[0];
							m_cosz[2]= m_cosx[0]*m_cosy[1] - m_cosx[1]*m_cosy[0];

						}
					}
				}

				leng=std::sqrt(m_cosx[2]*m_cosx[2] + m_cosy[2]*m_cosy[2] + m_cosz[2]*m_cosz[2]);

				m_cosx[2]= m_cosx[2]/leng;
				m_cosy[2]= m_cosy[2]/leng;
				m_cosz[2]= m_cosz[2]/leng;

				for (int i=0; i<3; i++){
					if (m_cosx[i] == 0.0) m_cosx[i] = 0.001;
					if (m_cosy[i] == 0.0) m_cosy[i] = 0.001;
					if (m_cosz[i] == 0.0) m_cosz[i] = 0.001;
				}


				//matrix inversion
				/*det = m_cosx[0]*( m_cosy[1]*m_cosz[2] - m_cosy[2]*m_cosz[1] ) - m_cosy[0]*( m_cosx[1]*m_cosz[2] - m_cosx[2]*m_cosz[1]) + 
				m_cosz[0]*( m_cosx[1]*m_cosy[2] - m_cosx[2]*m_cosy[1]);


				q[2] = (m_cosx[1]*m_cosy[2] - m_cosx[2]*m_cosy[1]) / det;
				q[1] = (m_cosx[2]*m_cosz[1] - m_cosx[1]*m_cosz[2]) / det;
				q[0] = (m_cosy[1]*m_cosz[2] - m_cosy[2]*m_cosz[1] ) / det;
				q[5] = (m_cosx[2]*m_cosy[0] - m_cosx[0]*m_cosy[2]) / det;
				q[4] = (m_cosx[0]*m_cosz[2] - m_cosx[2]*m_cosz[0]) / det;
				q[3] =  (m_cosy[2]*m_cosz[0] - m_cosy[0]*m_cosz[2] ) / det;
				q[8] = (m_cosx[0]*m_cosy[1] - m_cosx[1]*m_cosy[0]) / det;
				q[7] = (m_cosx[1]*m_cosz[0] - m_cosx[0]*m_cosz[1]) / det;
				q[6] = (m_cosy[0]*m_cosz[1] - m_cosy[1]*m_cosz[0] ) / det;*/


				//matrix inversion
				A = m_cosy[1] - m_cosx[1] * m_cosy[0] / m_cosx[0];
				//if(A==0) A=0.001;
				V = m_cosz[1] - m_cosx[1] * m_cosz[0] / m_cosx[0];
				C = m_cosz[2] - m_cosy[2] * V/A + m_cosx[2] * (V/A * m_cosy[0] / m_cosx[0] - m_cosz[0] / m_cosx[0]);
				//if(C==0) C=0.001;
				inversionMat[2] = (-m_cosx[2] / m_cosx[0] - m_cosx[2] / A * m_cosx[1] / m_cosx[0] * m_cosy[0] / m_cosx[0] + m_cosy[2] / A * m_cosx[1] / m_cosx[0]) / C;
				inversionMat[1] = -inversionMat[2] * V / A - m_cosx[1] / (m_cosx[0] * A);
				inversionMat[0] = 1.0 / m_cosx[0] - inversionMat[1] * m_cosy[0] / m_cosx[0] - inversionMat[2] * m_cosz[0] / m_cosx[0];
				inversionMat[5] = (-m_cosy[2] / A + m_cosx[2] / A * m_cosy[0] / m_cosx[0]) / C;
				inversionMat[4] = -inversionMat[5] * V / A + 1.0 / A;
				inversionMat[3] = -inversionMat[4] * m_cosy[0] / m_cosx[0] - inversionMat[5] * m_cosz[0] / m_cosx[0];
				inversionMat[8] = 1.0 / C;
				inversionMat[7] = -inversionMat[8] * V / A;
				inversionMat[6] = -inversionMat[7] * m_cosy[0] / m_cosx[0] - inversionMat[8] * m_cosz[0] / m_cosx[0];

			}


			RawImage ColorDeconvolution::separate_stains(const RawImage &source, double conv_matrix[9])
			{
				// initialize 3 output colour images
				sedeen::Size imageSize = source.size();
				std::vector<RawImage> binaryImages;

				for (int i=0; i<NumOfStains; i++){
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
					double Rlog = -((255.0*log((R+1)/255.0))/log255);
					double Glog = -((255.0*log((G+1)/255.0))/log255);
					double Blog = -((255.0*log((B+1)/255.0))/log255);

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


                        //Need a condition around this using m_applyThreshold
                        unsigned char binaryValue;
                        if (m_applyThreshold) {
                            binaryValue = ((Rlogscaled + Glogscaled + Blogscaled) > m_threshold) ? 255 : 0;
                        }
                        else {
                            binaryValue = ((Rlogscaled + Glogscaled + Blogscaled) > 0.0) ? 255 : 0;
                        }

						binaryImages.at(i).setValue(x, y, 0, binaryValue);
						binaryImages.at(i).setValue(x, y, 1, 0);
						binaryImages.at(i).setValue(x, y, 2, 0);
						binaryImages.at(i).setValue(x, y, 3, 255);
					}
				}

				if( m_DisplayOption == DisplayOptions::STAIN1 ){					
					return binaryImages[0];
				}

				if( m_DisplayOption == DisplayOptions::STAIN2 ){
					return binaryImages[1];
				}

				if( m_DisplayOption == DisplayOptions::STAIN3 ){
					return binaryImages[2];
				}

			}

			RawImage ColorDeconvolution::doProcessData(const RawImage &source)
			{
				//matrix inversion
				double  conv_matrix[9];
				computeMatrixInverse( conv_matrix );

				//Stain Separation and combination
				return separate_stains( source, conv_matrix);

			}


			const ColorSpace& ColorDeconvolution::doGetColorSpace() const
			{
				return m_colorSpace;
			}

			//bool ColorDeconvolution::saveToCSVfile(const std::string& fileName )
			//{
			//	std::ofstream file (fileName);
            //
			//	if( file.is_open())
			//	{
			//		file << "RegionOfInterest;"; 
			//		for (int i=0; i < 3; i++)
			//		{
			//			file << m_MODx[i] << ";" << m_MODy[i] <<
			//				";" << m_MODz[i] << ";"; 
			//		}
			//		file << std::endl;
			//	}
			//	else
			//		return 0;
			//	file.close();
			//	return 1;
			//}

			//void ColorDeconvolution::loadFromCSV( const std::string& fileName, const std::string& marker )
			//{
			//	std::ifstream file(fileName);
			//	std::vector<std::string> row;
			//	std::string line;
			//	std::string cell;

			//	int i=0;
			//	while( std::getline(file, line) )
			//	{
			//		std::stringstream lineStream(line);
			//		row.clear();
            //
			//		while( std::getline( lineStream, cell, ';' ) )
			//			row.push_back( cell );
            //
			//		int szdiv = (row.size()-1) % 3;
			//		if(!row.empty() && row[0] == marker && szdiv==0)
			//		{
			//			for(unsigned int j=1; j < row.size(); j+=3)
			//			{
			//				m_MODx[i] = stod(row.at(j));
			//				m_MODy[i] = stod(row.at(j+1));
			//				m_MODz[i]  = stod(row.at(j+2));
			//				i++;
			//			}
			//		}
			//	}
			//}

		} // namespace tile

		void getmeanRGBODfromROI(RawImage ROI, double rgbOD[3])
		{
			if( ROI.isNull() )
				return;

			int imageSize = ROI.size().width()*ROI.size().height();
			double log255= log(255.0);
			int y=0, x =0;
			for (int i=0; i < imageSize ; i++){
				x = i%ROI.size().width();
				if( x == 0 && i!=0)
					y++;

				// rescale to match original paper values
				rgbOD[0] =rgbOD[0] + (-((255.0*log(((ROI.at(x, y, 0).as<double>() +1)/255.0))/log255)));
				rgbOD[1] =rgbOD[1] + (-((255.0*log(((ROI.at(x, y, 1).as<double>() +1)/255.0))/log255)));
				rgbOD[2] =rgbOD[2] + (-((255.0*log(((ROI.at(x, y, 2).as<double>() +1)/255.0))/log255)));
			}
			rgbOD[0] = rgbOD[0] / imageSize;
			rgbOD[1] = rgbOD[1] / imageSize;
			rgbOD[2] = rgbOD[2] / imageSize;

		}//end getmeanRGBODfromROI

		void getStainsComponents(std::shared_ptr<tile::Factory> source,
			const std::vector<std::shared_ptr<GraphicItemBase>> region_of_interests,
			const Size &rescaled_resolutions, double conv_matrix[9])
		{

			for(int i=0; i < 9; i++)
				conv_matrix[i] = 0;

			// Get image from the output factory
			auto compositor = image::tile::Compositor(source);
			double  rgbOD[3];
			for(int j=0; j < 3; j++)
			{
				rgbOD[0]=0;
				rgbOD[1]=0;
				rgbOD[2]=0;
				Rect rect = containingRect(region_of_interests.at(j)->graphic());
				//RawImage ROI = compositor.getImage(rect, rescaled_resolutions);
				RawImage ROI = compositor.getImage(rect, Size(rect.width(), rect.height()));
				getmeanRGBODfromROI(ROI, rgbOD);

				conv_matrix[j*3] = rgbOD[0];
				conv_matrix[j*3+1] = rgbOD[1];
				conv_matrix[j*3+2] = rgbOD[2];
			}

		}//end getStainsComponents

	} // namespace image
} // namespace sedeen