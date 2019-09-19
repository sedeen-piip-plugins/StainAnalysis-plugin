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
#include "StainVectorFromSVD.h"

//OpenCV include
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

//Get processing time with C++11 chrono
#include <chrono>
#include <random>
#include <ctime>
#include <cmath>
#include <sstream>
#include <limits>

namespace sedeen {
namespace image {

    long long doSomethingWithSVD(std::shared_ptr<tile::Factory> source, const Size &rescaled_resolutions,
        double stainVec_matrix[9]) {
        //Initialize 64-bit random number generator
        std::random_device rd;
        std::mt19937_64 rgen(rd()); //64-bit Mersenne Twister


        //TEMP: don't modify it
        //Clear the stainVec_matrix
        //for (int i = 0; i < 9; i++) {
        //    stainVec_matrix[i] = 0;
        //}


        //This determines the sample size
        int numberOfSamplePixels = 1000;

        //Define OpenCV Mat structure with numberOfSamplePixels rows, RGB columns, elements are type double
        cv::Mat sampledPixelsMatrix(numberOfSamplePixels, 3, cv::DataType<double>::type);


        //temporary threshold on the average OD
        double avgODThreshold = 0.15;
        double percentileThreshold = 1.0;

        //Time the operation
        //Start the steady clock
        auto startTime = std::chrono::steady_clock::now();

        //Get features of the image
        s32 numResLevels = source->getNumLevels();
        //The highest resolution is always level 0
        s32 highResLevel = 0;
        s32 numHighResTiles = source->getNumTiles(highResLevel);
        s32 numTilePixels = static_cast<s32>(source->getTileSize().width() * source->getTileSize().height());
        //Some WSIs have multiple focus planes. Act on the default one
        auto numFocusPlanes = tile::getNumFocusPlanes(*source);
        auto defaultFocusPlane = tile::getDefaultFocusPlane(*source);
        //Act on the default band (Brightfield or Fluorescence)
        auto numBands = tile::getNumBands(*source);
        auto defaultBand = tile::getDefaultBand(*source);


        //Create an initialized array to store the number of required pixels from each tile
        std::unique_ptr<u16[]> tileSamplingCountArray = std::make_unique<u16[]>(numHighResTiles);
        //Initialize a random distribution to choose tile indices
        std::uniform_int_distribution<s32> randTileIndex(0, numHighResTiles - 1);
        for (int spx = 0; spx < numberOfSamplePixels; spx++) {
            tileSamplingCountArray[randTileIndex(rgen)]++;
        }


#include <fstream>
        //Temp file output
        std::fstream tempOut;
        tempOut.open("D:\\mschumaker\\projects\\Sedeen\\testData\\output\\tempout.txt", std::fstream::out);


        //Create a TileServer to be able to access tiles from the factory
        std::unique_ptr<tile::TileServer> theTileServer
            = std::unique_ptr<tile::TileServer>(new tile::TileServer(source));

        //Loop over the tiles in the high res image
        long numPixelsAddedToMatrix = 0;
        for (int tl = 0; tl < numHighResTiles; tl++) {
            if (tileSamplingCountArray[tl] > 0) {




                //Create an array of pixel indices
                std::unique_ptr<u8[]> pixelSamplingArray = std::make_unique<u8[]>(numTilePixels);
                //Initialize a uniform random distribution
                std::uniform_int_distribution<s32> randPixelIndex(0, numTilePixels - 1);
                //Fill the array with the number of required pixels, no duplication
                for (int tpx = 0; tpx < tileSamplingCountArray[tl]; tpx++) {
                    int countLimit = 2 * numTilePixels; //Kind of high, but shouldn't be needed
                    bool freeLocationFound = false;
                    int attemptNumber = 0;
                    while (!freeLocationFound && (attemptNumber < countLimit)) {
                        s32 newPixelIndex = randPixelIndex(rgen);
                        if (pixelSamplingArray[newPixelIndex] == 0) {
                            pixelSamplingArray[newPixelIndex] = 1;
                            freeLocationFound = true;
                        }
                        else if (pixelSamplingArray[newPixelIndex] == 1) {
                            freeLocationFound = false;
                            attemptNumber++;
                        }
                        else {
                            //invalid value. TODO: error handling
                            freeLocationFound = false;
                            attemptNumber++;
                        }
                    }
                }

                //Retrieve this tile, place in a RawImage so that pixel values are accessible
                auto tileIndex = tile::getTileIndex(*source, highResLevel, tl, defaultFocusPlane, defaultBand);
                RawImage tileImage = theTileServer->getTile(tileIndex);
                auto numPixels = tileImage.width() * tileImage.height();
                auto numChannels = sedeen::image::channels(tileImage);
                auto numElements = numPixels * numChannels;
                //Get the pixel order of the image: Interleaved or Planar
                PixelOrder pixelOrder = tileImage.order();

                //For every chosen pixel set to 1 in pixelSamplingArray
                for (int px = 0; px < numPixels; px++) {
                    if (pixelSamplingArray[px] == 1) {
                        double rgbOD[3] = { 0.0 };
                        unsigned int Rindex, Gindex, Bindex;
                        if (pixelOrder == PixelOrder::Interleaved) {
                            //RGB RGB RGB ... (if numChannels=3)
                            Rindex = px * numChannels + 0;
                            Gindex = px * numChannels + 1;
                            Bindex = px * numChannels + 2;
                        }
                        else if (pixelOrder == PixelOrder::Planar) {
                            //RRR... GGG... BBB...
                            Rindex = 0 * numPixels + px;
                            Gindex = 1 * numPixels + px;
                            Bindex = 2 * numPixels + px;
                        }
                        else {
                            //Invalid value of pixelOrder
                            break;
                        }
                        //Check the values
                        if ((Rindex >= numElements) || (Gindex >= numElements) || (Bindex >= numElements)) {
                            break;
                        }
                        //Get the optical density values
                        rgbOD[0] = StainVectorMath::ConvertRGBtoOD(static_cast<double>((tileImage[Rindex]).as<s32>()));
                        rgbOD[1] = StainVectorMath::ConvertRGBtoOD(static_cast<double>((tileImage[Gindex]).as<s32>()));
                        rgbOD[2] = StainVectorMath::ConvertRGBtoOD(static_cast<double>((tileImage[Bindex]).as<s32>()));

                        if (rgbOD[0] + rgbOD[1] + rgbOD[2] > avgODThreshold) {
                            sampledPixelsMatrix.at<double>(numPixelsAddedToMatrix, 0) = rgbOD[0];
                            sampledPixelsMatrix.at<double>(numPixelsAddedToMatrix, 1) = rgbOD[1];
                            sampledPixelsMatrix.at<double>(numPixelsAddedToMatrix, 2) = rgbOD[2];
                            numPixelsAddedToMatrix++;

                            //Let's write these to a file to check that something's happening
                            //std::stringstream sv;
                            //sv << numPixelsAddedToMatrix-1 << ") Tile " << tl << ", Pixel " << px << ": ";
                            //sv << sampledPixelsMatrix.at<double>(numPixelsAddedToMatrix-1, 0) << ", ";
                            //sv << sampledPixelsMatrix.at<double>(numPixelsAddedToMatrix-1, 1) << ", ";
                            //sv << sampledPixelsMatrix.at<double>(numPixelsAddedToMatrix-1, 2);
                            //tempOut << sv.str() << std::endl;
                        }
                    }
                }

                //Free the pixelSamplingArray
                pixelSamplingArray.release();
                //tileImage should go out of scope
            }
        }//end for each tile


        //Time to get the covariance matrix
        //Shorten the sampledPixelsMatrix
        sampledPixelsMatrix.resize(numPixelsAddedToMatrix);


        //Get column means. Use those as the mean input to calcCovarMatrix


        //use PCA to simplify math
        cv::PCA pcAnalysis(sampledPixelsMatrix, cv::noArray(), cv::PCA::DATA_AS_ROW, 2); //2 components to project onto a plane
        //Project all of the sampled pixels into the new basis
        cv::Mat projectedPoints;
        pcAnalysis.project(sampledPixelsMatrix, projectedPoints);

        std::stringstream scov;
        scov << "The PCA eigenvalues are: " << pcAnalysis.eigenvalues << std::endl;
        scov << "The PCA basis vectors are: " << pcAnalysis.eigenvectors << std::endl;
        scov << "The first projected point is: " << projectedPoints.row(0) << std::endl;
        //scov << "What's different I think is the mean array that's used: " << pcAnalysis.mean << std::endl;


        //Is projecting into the new basis all that difficult?
        cv::Mat basisVectors;
        cv::transpose(pcAnalysis.eigenvectors, basisVectors);


        if (basisVectors.at<double>(0, 0) < 0.0) {
            basisVectors.at<double>(0, 0) *= -1.0;
            basisVectors.at<double>(1, 0) *= -1.0;
            basisVectors.at<double>(2, 0) *= -1.0;
        }
        if (basisVectors.at<double>(0, 1) < 0.0) {
            basisVectors.at<double>(0, 1) *= -1.0;
            basisVectors.at<double>(1, 1) *= -1.0;
            basisVectors.at<double>(2, 1) *= -1.0;
        }

        scov << "The basisVectors: " << std::endl;
        scov << basisVectors << std::endl;


        cv::Mat projectedPointsByMM;
        cv::gemm(sampledPixelsMatrix, basisVectors, 1, cv::Mat(), 0, projectedPointsByMM, 0);

        scov << "The gemm output: " << std::endl;
        scov << projectedPointsByMM << std::endl;

        //Get angular coordinates with respect to the eigenvectors (the new orthogonal basis)
        cv::Mat angularCoords(projectedPointsByMM.rows, 1, cv::DataType<float>::type);

        for (auto p = angularCoords.begin<float>(); p != angularCoords.end<float>(); p++) {
            int row = static_cast<int>(p.lpos());
            bool angleUndef = (projectedPointsByMM.at<float>(row, 0) == 0.0)
                && (projectedPointsByMM.at<float>(row, 1) == 0.0);
            //Use maximum float value as undefined value
            float undefined = std::numeric_limits<float>::max();
            *p = angleUndef ? undefined
                : std::atan2(projectedPointsByMM.at<float>(row, 1), projectedPointsByMM.at<float>(row, 0));
        }

        //Create a histogram of values between -pi and pi
        cv::Mat angleHist;
        int channels[1] = { 0 };
        int histSize[1] = { 1024 };
        float pi = static_cast<float>(CV_PI);
        float range[] = { -pi, pi };
        const float* histRange[] = { range };
        bool uniform = true;
        bool accumulate = false;
        cv::calcHist(&angularCoords, 1, channels, cv::Mat(), angleHist, 1, histSize, histRange, uniform, accumulate);
        //The resulting type of angleHist elements is float

        scov << "This is angleHist size: " << std::endl;
        scov << (angleHist.size)[0] << std::endl;
        scov << angleHist.type() << std::endl;

        //Count how many elements were added to the histogram
        cv::Scalar scalarTotal = cv::sum(angleHist);
        double histoCountTotal = scalarTotal[0];

        //Calculate the slope and intercept to convert from bin to angle
        //angle = intercept + slope*bin
        double intercept = static_cast<double>(range[0]);
        double slope = static_cast<double>(range[1] - range[0]) / static_cast<double>(histSize[0]);

        //Find the bins containing the lower and upper percentiles
        double lowerFraction = (percentileThreshold) / 100.0;
        double upperFraction = (100.0 - percentileThreshold) / 100.0;
        bool lowerPassed(false), upperPassed(false);
        double lowerBin(-1.0), upperBin(-1.0);
        double cumulativeSum(0.0);
        for (int bin = 0; bin < (angleHist.size)[0]; bin++) {
            double prevFraction = cumulativeSum / histoCountTotal;
            cumulativeSum += static_cast<double>(angleHist.at<float>(bin, 0));
            double currentFraction = cumulativeSum / histoCountTotal;

            //scov << bin << " currentFraction: " << currentFraction << std::endl;

            if (!lowerPassed && (currentFraction >= lowerFraction)) {
                lowerPassed = true;
                //linearly interpolate a bin position
                lowerBin = bin + (lowerFraction - prevFraction) / (currentFraction - prevFraction);
            }
            if (!upperPassed && (currentFraction >= upperFraction)) {
                upperPassed = true;
                //linearly interpolate a bin position
                upperBin = bin + (upperFraction - prevFraction) / (currentFraction - prevFraction);
            }
            if (lowerPassed && upperPassed) { break; }
        }
        //Convert from bin values to angles
        double lowerAngle = intercept + slope * lowerBin;
        double upperAngle = intercept + slope * upperBin;


        scov << "Lower bin number: " << lowerBin << ", upper bin number: " << upperBin << std::endl;
        scov << "Lower angle: " << lowerAngle << ", upper angle: " << upperAngle << std::endl;

        //Convert to Cartesian on the PCA projected plane
        cv::Mat stainVectorProjections(2, 2, cv::DataType<double>::type);
        stainVectorProjections.at<double>(0, 0) = std::cos(lowerAngle);
        stainVectorProjections.at<double>(0, 1) = std::sin(lowerAngle);
        stainVectorProjections.at<double>(1, 0) = std::cos(upperAngle);
        stainVectorProjections.at<double>(1, 1) = std::sin(upperAngle);

        //Back-project to get un-normalized stain vectors
        cv::Mat stainVectorOutput;
        pcAnalysis.backProject(stainVectorProjections, stainVectorOutput);



        //pick up here!!!  Do the backprojection using my basisVectors




        
        scov << "Well, here are what it calculates as the stain vectors:" << std::endl;
        scov << stainVectorOutput;




        tempOut << scov.str() << std::endl;
        tempOut.close();




        //How long did it take to get the basis vectors?
        auto afterBasisVectorsTime = std::chrono::steady_clock::now();

        //How long did getting the basis vectors take?
        long long gettingBasisVectorsDuration = std::chrono::duration_cast<std::chrono::microseconds>(afterBasisVectorsTime - startTime).count();





        //Prep for returning
        tileSamplingCountArray.release();



        //return numHighResTiles;

        return gettingBasisVectorsDuration;




        //return numPixelsAddedToMatrix;
    }//end SVD



} // namespace image
} // namespace sedeen
