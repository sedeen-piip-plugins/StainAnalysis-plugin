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

#ifndef SEDEEN_SRC_FILTER_STAINVECTORNIETHAMMER_H
#define SEDEEN_SRC_FILTER_STAINVECTORNIETHAMMER_H

#include "Global.h"
#include "Geometry.h"
#include "Image.h"

#include "StainVectorOpenCV.h"

namespace sedeen {
namespace image {

class PATHCORE_IMAGE_API StainVectorNiethammer : public StainVectorOpenCV {
public:
    StainVectorNiethammer(std::shared_ptr<tile::Factory> source);
    ~StainVectorNiethammer();

    ///Fill the 9-element array with three stain vectors
    virtual void ComputeStainVectors(double outputVectors[9]);
    ///Overload of the basic method, includes parameters needed by the algorithm
    void ComputeStainVectors(double outputVectors[9], const int sampleSize, 
        const double ODthreshold = 0.15, const double percentileThreshold = 1.0);

    ///Get/Set the average optical density threshold
    inline const double GetODThreshold() const { return m_avgODThreshold; }
    ///Get/Set the average optical density threshold
    inline void SetODThreshold(const double t) { m_avgODThreshold = t; }

    ///Get/Set the percentile threshold
    inline const double GetPercentileThreshold() const { return m_percentileThreshold; }
    ///Get/Set the percentile threshold
    inline void SetPercentileThreshold(const double p) { m_percentileThreshold = p; }

    ///Get/Set the sample size, the number of pixels to choose
    inline const int GetSampleSize() const { return m_sampleSize; }
    ///Get/Set the sample size, the number of pixels to choose
    inline void SetSampleSize(const int s) { m_sampleSize = s; }

private:
    double m_avgODThreshold;
    double m_percentileThreshold;

    ///The number of pixels that should be used to calculate the stain vectors
    int m_sampleSize;
};

} // namespace image
} // namespace sedeen
#endif
