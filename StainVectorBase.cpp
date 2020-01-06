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

#include "StainVectorBase.h"

namespace sedeen {
namespace image {

StainVectorBase::StainVectorBase(std::shared_ptr<tile::Factory> source) 
    : m_sourceFactory(source),
    m_randomWSISampler(std::make_shared<RandomWSISampler>(source))
{
}//end constructor

StainVectorBase::~StainVectorBase(void) {
}//end destructor

void StainVectorBase::ComputeStainVectors(double (&outputVectors)[9]) {
}//end ComputeStainVectors

bool StainVectorBase::AreEqual(cv::InputArray array1, cv::InputArray array2) {
    // treat two empty arrays as identical
    if (array1.empty() && array2.empty()) {
        return true;
    }
    // if dimensionality is not identical, these arrays are not identical
    if (array1.cols() != array2.cols() || array1.rows() != array2.rows() || array1.dims() != array2.dims()) {
        return false;
    }
    //Compare NOT equal, then count NON-zero (there isn't a countZero function in OpenCV).
    cv::Mat diff;
    cv::compare(array1, array2, diff, cv::CmpTypes::CMP_NE);
    int nz = cv::countNonZero(diff);
    return (nz == 0);
}//end AreEqual

} // namespace image
} // namespace sedeen
