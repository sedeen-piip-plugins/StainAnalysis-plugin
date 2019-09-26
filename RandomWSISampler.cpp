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

#include "RandomWSISampler.h"

//For now, include StainVectorMath here, but try to do the OD conversion and thresholding
//in a kernel, and use a factory to apply it before passing the factory to this class
#include "StainVectorMath.h"

#include <fstream>
#include <sstream>

namespace sedeen {
namespace image {

RandomWSISampler::RandomWSISampler(std::shared_ptr<tile::Factory> source) 
    : m_sourceFactory(source),
    m_rgen((std::random_device())()) //Initialize random number generation
{
}//end constructor

RandomWSISampler::~RandomWSISampler(void) {
}//end destructor


bool RandomWSISampler::ChooseRandomPixels(cv::Mat outputMatrix, int numberOfPixels, double ODthreshold,
    int level /* = 0 */, int focusPlane /* = -1 */, int band /* = -1 */) {
    assert(nullptr != GetSourceFactory());
    auto source = this->GetSourceFactory();
    //get info about the whole slide image from the source factory
    s32 numResolutionLevels = source->getNumLevels();
    //Some WSIs have multiple focal planes. The WSI will have a default focus plane set
    auto numFocusPlanes = tile::getNumFocusPlanes(*source);
    auto defaultFocusPlane = tile::getDefaultFocusPlane(*source);
    //Bands (Brightfield or Fluorescence)
    auto numBands = tile::getNumBands(*source);
    auto defaultBand = tile::getDefaultBand(*source);

    //Check the level, focusPlane, and band argument values
    //Highest resolution level is 0, level must be within range
    if ((level < 0) || (level >= numResolutionLevels)) { return false; }
    //If focusPlane is -1, choose the default focus plane
    s32 chosenFocusPlane = static_cast<s32>((focusPlane < 0) ? defaultFocusPlane : focusPlane);
    if (chosenFocusPlane >= numFocusPlanes) { return false; }
    //If band is -1, choose the default band
    s32 chosenBand = static_cast<s32>((band < 0) ? defaultBand : band);
    if (chosenBand >= numBands) { return false; }

    //Get the number of tiles on the chosen level, and the pixels per tile
    s32 numTilesOnLevel = source->getNumTiles(level);
    s32 numTilePixels = static_cast<s32>(source->getTileSize().width() * source->getTileSize().height());


    //Temp file output
    std::fstream tempOut;
    tempOut.open("D:\\mschumaker\\projects\\Sedeen\\testData\\output\\tempout.txt", std::fstream::out);

    tempOut << "We're at least getting into ChooseRandomPixels" << std::endl;




    //Create an initialized array to store the number of required pixels from each tile
    std::unique_ptr<u16[]> tileSamplingCountArray = std::make_unique<u16[]>(numTilesOnLevel);
    //Initialize a uniform distribution to choose tile indices
    std::uniform_int_distribution<s32> randTileIndex(0, numTilesOnLevel - 1);
    for (int spx = 0; spx < numberOfPixels; spx++) {
        s32 newIndex = randTileIndex(m_rgen);
        tileSamplingCountArray[newIndex]++;
        tempOut << newIndex << std::endl;
    }

    //Create a TileServer to be able to access tiles from the factory
    std::unique_ptr<tile::TileServer> theTileServer
        = std::unique_ptr<tile::TileServer>(new tile::TileServer(source));

    //Define OpenCV Mat structure with numberOfSamplePixels rows, RGB columns, elements are type double
    cv::Mat sampledPixelsMatrix(numberOfPixels, 3, cv::DataType<double>::type);








    //Loop over the tiles in the high res image
    long numPixelsAddedToMatrix = 0;
    for (int tl = 0; tl < numTilesOnLevel; tl++) {
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
                    s32 newPixelIndex = randPixelIndex(m_rgen);
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
            auto tileIndex = tile::getTileIndex(*source, level, tl, defaultFocusPlane, defaultBand);
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

                    if (rgbOD[0] + rgbOD[1] + rgbOD[2] > ODthreshold) {
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





    tempOut.close();











    //Copy the pixel data from the temp matrix to the output matrix
    outputMatrix = cv::Mat(numPixelsAddedToMatrix, 3, cv::DataType<double>::type);




    return true;
}//end ChooseRandomPixels



} // namespace image
} // namespace sedeen
