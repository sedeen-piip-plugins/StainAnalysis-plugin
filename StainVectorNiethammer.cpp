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

#include "StainVectorNiethammer.h"

#include <chrono>
#include <random>

 //For now, include StainVectorMath here, but try to do the OD conversion and thresholding
 //in a kernel, and use a factory to apply it before passing the factory to this class
#include "StainVectorMath.h"

namespace sedeen {
namespace image {

StainVectorNiethammer::StainVectorNiethammer(std::shared_ptr<tile::Factory> source)
    : StainVectorBase(source), m_sourceFactory(source), m_avgODThreshold(0.15)
{
    //Initialize 64-bit random number generator
    std::random_device rd;
    std::mt19937_64 rgen(rd()); //64-bit Mersenne Twister

}//end constructor

StainVectorNiethammer::~StainVectorNiethammer(void) {
}//end destructor


long long StainVectorNiethammer::ChooseRandomPixels(cv::Mat outputMatrix, int numberOfPixels, bool suppressZeros) {
    assert(nullptr != m_sourceFactory);


    return 0;
}//end ChooseRandomPixels



} // namespace image
} // namespace sedeen
