#include "ColorDeconvolutionKernel.h"

namespace sedeen {
	namespace image {
		namespace tile {
			ColorDeconvolution::ColorDeconvolution( Behavior behavior, DisplyOptions displyOption, double conv_matrix[9], double threshold, const std::string& path_to_root):
				m_threshold(threshold),
				count(0),
				m_Stain(behavior),
				log255(std::log(255.0)),
				m_displyOption(displyOption),
        m_colorSpace(ColorModel::RGBA, ChannelType::UInt8)

			{
				SetStainMatrice(m_Stain, conv_matrix, path_to_root);
				setDisplayOptions(m_displyOption);

				/*std::string name = "C:/sedeen/logfile1.txt";
				log_file= std::ofstream(name, std::ios_base::out | std::ios_base::app );*/
			}


			ColorDeconvolution::~ColorDeconvolution(void)
			{
			}

			void ColorDeconvolution::setDisplayOptions(DisplyOptions displyOption)
			{
				if( displyOption == DisplyOptions::STAIN1 )
					m_displyOption = DisplyOptions::STAIN1;

				if( displyOption == DisplyOptions::STAIN2 )
					m_displyOption = DisplyOptions::STAIN2;

				if( displyOption == DisplyOptions::STAIN3 )
					m_displyOption = DisplyOptions::STAIN3;

			}


			void ColorDeconvolution::SetStainMatrice(Behavior behavior, double conv_matrix[9], const std::string& path_to_root)
			{
				path dir_path = path_to_root +"sedeen/";
				if(behavior == Behavior::RegionOfInterest)
				{
					for (int i=0; i < 3; i++)
					{
						MODx[i]= conv_matrix[i*3];
						MODy[i]= conv_matrix[i*3+1];
						MODz[i]= conv_matrix[i*3+2];
					}

					path dir_path = path_to_root +"sedeen/";
					std::string file_name = "StainsFile.csv";

					/*path dir_path = "C:/sedeen/";
					std::string file_name = "StainsFile.csv";*/
					try
					{
						//&& boost::filesystem::is_directory(dir_path)
						if( !boost::filesystem::exists( dir_path ) )
						{
							if(boost::filesystem::create_directory( dir_path ))
							{
								std::string fn = dir_path.string() + file_name;
								saveToCSVfile(fn);
							}

						}
						else
						{
							std::string fn = dir_path.string() + file_name;
							saveToCSVfile(fn);
						}

					}
					catch(boost::filesystem::filesystem_error const & e)
					{

					}

				}

				if(behavior == Behavior::LoadFromFile)
				{
					std::string file_name = path_to_root +"sedeen/StainsFile.csv";
					//std::string fileName = "C:/sedeen/StainsFile.csv";
					loadFromCSV(file_name, "RegionOfInterest");
				}

				if (behavior == Behavior::HematoxylinPEosin){
					//path dir_path = "C:/sedeen/";
					path dir_path = path_to_root +"sedeen/";
					std::string file_name = "DefaultStainsFile.csv";
					path path_found = path();

					bool thereIsAny = image::serachForfile( dir_path, file_name, path_found );
					if(thereIsAny)
					{
						std::string fn = dir_path.string() + file_name;
						loadFromCSV(fn, "HematoxylinPEosin");
					}
					else
					{
						// GL Haem matrix
						MODx[0]= 0.644211; //0.650;
						MODy[0]= 0.716556; //0.704;
						MODz[0]= 0.266844; //0.286;
						// GL Eos matrix
						MODx[1]= 0.092789; //0.072;
						MODy[1]= 0.954111; //0.990;
						MODz[1]= 0.283111; //0.105;
						// Zero matrix
						MODx[2]= 0.0;
						MODy[2]= 0.0;
						MODz[2]= 0.0;
					}

				}

				if (behavior == Behavior::HematoxylinPDAB){

					//path dir_path = "C:/sedeen/";
					path dir_path = path_to_root +"sedeen/";
					std::string file_name = "DefaultStainsFile.csv";
					path path_found = path();

					bool thereIsAny = image::serachForfile( dir_path, file_name, path_found );
					if(thereIsAny)
					{
						std::string fn = dir_path.string() + file_name;
						loadFromCSV(fn, "HematoxylinPDAB");
					}
					else
					{
						// 3,3-diamino-benzidine tetrahydrochloride
						// Haem matrix
						MODx[0]= 0.650;
						MODy[0]= 0.704;
						MODz[0]= 0.286;
						// DAB matrix
						MODx[1]= 0.268;
						MODy[1]= 0.570;
						MODz[1]= 0.776;
						// Zero matrix
						MODx[2]= 0.0;
						MODy[2]= 0.0;
						MODz[2]= 0.0;
					}
				}

				if (behavior == Behavior::HematoxylinPEosinPDAB){
					//path dir_path = "C:/sedeen/";
					path dir_path = path_to_root +"sedeen/";
					std::string file_name = "DefaultStainsFile.csv";
					path path_found = path();

					bool thereIsAny = image::serachForfile( dir_path, file_name, path_found );
					if(thereIsAny)
					{
						std::string fn = dir_path.string() + file_name;
						loadFromCSV(fn, "HematoxylinPEosinPDAB");
					}
					else
					{
						// Haem matrix
						MODx[0]= 0.650;
						MODy[0]= 0.704;
						MODz[0]= 0.286;
						// Eos matrix
						MODx[1]= 0.072;
						MODy[1]= 0.990;
						MODz[1]= 0.105;
						// DAB matrix
						MODx[2]= 0.268;
						MODy[2]= 0.570;
						MODz[2]= 0.776;
					}
				}

			}

			void ColorDeconvolution::computeMatrixInvers( double inversionMat[9] )
			{
				double leng, A, V, C;
				double  len[NumOfStains];
				//double  q[9];

				for (int i=0; i<NumOfStains; i++){
					//normalise vector length
					cosx[i]=cosy[i]=cosz[i]=0.0;
					len[i]=std::sqrt(MODx[i]*MODx[i] + MODy[i]*MODy[i] + MODz[i]*MODz[i]);
					if (len[i] != 0.0){
						cosx[i]= MODx[i]/len[i];
						cosy[i]= MODy[i]/len[i];
						cosz[i]= MODz[i]/len[i];
					}
				}

				// translation matrix
				if (cosx[1]==0.0){ //2nd colour is unspecified
					if (cosy[1]==0.0){
						if (cosz[1]==0.0){
							cosx[1]=cosz[0];
							cosy[1]=cosx[0];
							cosz[1]=cosy[0];
						}
					}
				}


				if (cosx[2]==0.0){ // 3rd colour is unspecified
					if (cosy[2]==0.0){
						if (cosz[2]==0.0){

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

				for (int i=0; i<3; i++){
					if (cosx[i] == 0.0) cosx[i] = 0.001;
					if (cosy[i] == 0.0) cosy[i] = 0.001;
					if (cosz[i] == 0.0) cosz[i] = 0.001;
				}


				//matrix inversion
				/*det = cosx[0]*( cosy[1]*cosz[2] - cosy[2]*cosz[1] ) - cosy[0]*( cosx[1]*cosz[2] - cosx[2]*cosz[1]) + 
				cosz[0]*( cosx[1]*cosy[2] - cosx[2]*cosy[1]);


				q[2] = (cosx[1]*cosy[2] - cosx[2]*cosy[1]) / det;
				q[1] = (cosx[2]*cosz[1] - cosx[1]*cosz[2]) / det;
				q[0] = (cosy[1]*cosz[2] - cosy[2]*cosz[1] ) / det;
				q[5] = (cosx[2]*cosy[0] - cosx[0]*cosy[2]) / det;
				q[4] = (cosx[0]*cosz[2] - cosx[2]*cosz[0]) / det;
				q[3] =  (cosy[2]*cosz[0] - cosy[0]*cosz[2] ) / det;
				q[8] = (cosx[0]*cosy[1] - cosx[1]*cosy[0]) / det;
				q[7] = (cosx[1]*cosz[0] - cosx[0]*cosz[1]) / det;
				q[6] = (cosy[0]*cosz[1] - cosy[1]*cosz[0] ) / det;*/


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

						unsigned char binaryValue = ( (Rlogscaled + Glogscaled + Blogscaled)  > m_threshold) ? 255:0;						
						binaryImages.at(i).setValue(x, y, 0, binaryValue);
						binaryImages.at(i).setValue(x, y, 1, 0);
						binaryImages.at(i).setValue(x, y, 2, 0);
						binaryImages.at(i).setValue(x, y, 3, 255);
					}
				}

				if( m_displyOption == DisplyOptions::STAIN1 ){					
					return binaryImages[0];
				}

				if( m_displyOption == DisplyOptions::STAIN2 ){
					return binaryImages[1];
				}

				if( m_displyOption == DisplyOptions::STAIN3 ){
					return binaryImages[2];
				}

			}

			RawImage ColorDeconvolution::doProcessData(const RawImage &source)
			{
				//matrix inversion
				double  conv_matrix[9];
				computeMatrixInvers( conv_matrix );

				//Stain Separation and combination
				return separate_stains( source, conv_matrix);

			}


			const ColorSpace& ColorDeconvolution::doGetColorSpace() const
			{
				return m_colorSpace;
			}

			bool ColorDeconvolution::saveToCSVfile(const std::string& fileName )
			{
				std::ofstream file (fileName);

				if( file.is_open())
				{
					file << "RegionOfInterest;"; 
					for (int i=0; i < 3; i++)
					{
						file << MODx[i] << ";" << MODy[i] <<
							";" << MODz[i] << ";"; 
					}
					file << std::endl;
				}
				else
					return 0;
				file.close();
				return 1;
			}

			void ColorDeconvolution::loadFromCSV( const std::string& fileName, const std::string& marker )
			{
				std::ifstream file(fileName);
				std::vector<std::string> row;
				std::string line;
				std::string cell;

				int i=0;
				while( std::getline(file, line) )
				{
					std::stringstream lineStream(line);
					row.clear();

					while( std::getline( lineStream, cell, ';' ) )
						row.push_back( cell );

					int szdiv = (row.size()-1) % 3;
					if(!row.empty() && row[0] == marker && szdiv==0)
					{
						for(unsigned int j=1; j < row.size(); j+=3)
						{
							MODx[i] = stod(row.at(j));
							MODy[i] = stod(row.at(j+1));
							MODz[i]  = stod(row.at(j+2));
							i++;
						}
					}
				}
			}

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

		}

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

		}

		bool serachForfile( const path & dir_path, const std::string & file_name, path & path_found )
		{
			if ( !exists( dir_path ) ) 
				return false;

			directory_iterator end_itr; // default construction yields past-the-end

			for ( directory_iterator itr( dir_path ); itr != end_itr; ++itr )
			{
				if ( is_directory(itr->status()) )
				{
					if ( serachForfile( itr->path(), file_name, path_found ) ) 
						return true;
				}
				else if ( itr->path().filename() == file_name )
				{
					path_found = itr->path();
					return true;
				}
			}
			return false;
		}

		void setImagePath( const std::string& path )
		{

		}


	} // namespace image
} // namespace sedeen