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

#ifndef SEDEEN_SRC_FILTER_RANDOMWSISAMPLER_H
#define SEDEEN_SRC_FILTER_RANDOMWSISAMPLER_H

#include "Global.h"
#include "Geometry.h"
#include "Image.h"

#include <chrono>
#include <random>

//OpenCV include
#include <opencv2/core/core.hpp>

namespace sedeen {
namespace image {

class PATHCORE_IMAGE_API RandomWSISampler {

public:
    RandomWSISampler(std::shared_ptr<tile::Factory> source);
    virtual ~RandomWSISampler();

    ///Populate an OutputArray with pixels chosen without duplication from the source tile factory
    virtual bool ChooseRandomPixels(cv::OutputArray outputMatrix, const long int numberOfPixels, const double ODthreshold,
        const int level = 0, const int focusPlane = -1, const int band = -1); //Negative indicates to use the source default values

    //TODO: make a random tile and sub-tile chooser too

protected:
    ///Allow derived classes to get the source factory pointer
    inline std::shared_ptr<tile::Factory> GetSourceFactory() { return m_sourceFactory; }
    ///Allow derived classes access to the random number generator (64-bit Mersenne Twister)
    std::mt19937_64 m_rgen; //64-bit Mersenne Twister

private:
    std::shared_ptr<tile::Factory> m_sourceFactory;

};

} // namespace image
} // namespace sedeen
#endif
