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

#include "StainVectorMacenko.h"

#include <chrono>
#include <random>

 //For now, include StainVectorMath here, but try to do the OD conversion and thresholding
 //in a kernel, and use a factory to apply it before passing the factory to this class
#include "StainVectorMath.h"

namespace sedeen {
namespace image {

StainVectorMacenko::StainVectorMacenko(std::shared_ptr<tile::Factory> source)
    : StainVectorBase(source)
{


}//end constructor

StainVectorMacenko::~StainVectorMacenko(void) {
}//end destructor

void StainVectorMacenko::ComputeStainVectors(double outputVectors[9]) {

}


void StainVectorMacenko::ComputeStainVectors(double outputVectors[9], int sampleSize,
    double ODthreshold /* = 0.15 */, double percentileThreshold /* = 1.0 */) {
    assert(nullptr != this->GetSourceFactory());


    //Sample a set of pixel values from the source
    cv::Mat samplePixels;
    auto theSampler = this->GetRandomWSISampler();
    assert(nullptr != theSampler);
    bool samplingSuccess = theSampler->ChooseRandomPixels(samplePixels, sampleSize, ODthreshold);
    //move this to the other overload

}//end ComputeStainVectors




} // namespace image
} // namespace sedeen
