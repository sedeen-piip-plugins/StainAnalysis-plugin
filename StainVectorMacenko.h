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

#ifndef SEDEEN_SRC_FILTER_STAINVECTORMACENKO_H
#define SEDEEN_SRC_FILTER_STAINVECTORMACENKO_H

#include "Global.h"
#include "Geometry.h"
#include "Image.h"

#include "StainVectorBase.h"

 //OpenCV include
#include <opencv2/core/core.hpp>

namespace sedeen {
namespace image {

class PATHCORE_IMAGE_API StainVectorMacenko : StainVectorBase {
public:
    StainVectorMacenko(std::shared_ptr<tile::Factory> source);
    ~StainVectorMacenko();

    long long ChooseRandomPixels(cv::Mat outputMatrix, int numberOfPixels, bool suppressZeros);

private:
    std::shared_ptr<tile::Factory> m_sourceFactory;

private:
    //TODO: replace this with a threshold factory applied BEFORE creating an instance of this class
    double m_avgODThreshold;
};

} // namespace image
} // namespace sedeen
#endif
