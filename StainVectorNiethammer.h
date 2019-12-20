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

#include <array>

#include "Global.h"
#include "Geometry.h"
#include "Image.h"

#include "StainVectorOpenCV.h"

namespace sedeen {
namespace image {

class PATHCORE_IMAGE_API StainVectorNiethammer : public StainVectorOpenCV {
public:
    StainVectorNiethammer(std::shared_ptr<tile::Factory> source, double ODthreshold = 0.15, 
        double percentileThreshold = 1.0, double qAdjustmentFactor = 0.15);
    virtual ~StainVectorNiethammer();

    ///Fill the 9-element array with three stain vectors
    virtual void ComputeStainVectors(double (&outputVectors)[9]);
    ///Overload of the basic method, includes parameters needed by the algorithm
    void ComputeStainVectors(double (&outputVectors)[9], const double (&inputPriors)[9], const int sampleSize);
    ///Additional overload that does not require an array of prior stain vector values
    void ComputeStainVectors(double (&outputVectors)[9], const int sampleSize);

    ///Get/Set the average optical density threshold
    inline const double GetODThreshold() const { return m_avgODThreshold; }
    ///Get/Set the average optical density threshold
    inline void SetODThreshold(const double t) { m_avgODThreshold = t; }

    ///Get/Set the percentile threshold
    inline const double GetPercentileThreshold() const { return m_percentileThreshold; }
    ///Get/Set the percentile threshold
    inline void SetPercentileThreshold(const double pt) { m_percentileThreshold = pt; }

    ///Get/Set the q vector mixing ratio between s1 and s2 to get q1, q2
    inline const double GetQVectorMixRatio() const { return m_qVectorMixRatio; }
    ///Get/Set the q vector mixing ratio between s1 and s2 to get q1, q2
    inline void SetQVectorMixRatio(const double q) { m_qVectorMixRatio = q; }

    ///Get/Set the sample size, the number of pixels to choose
    inline const int GetSampleSize() const { return m_sampleSize; }
    ///Get/Set the sample size, the number of pixels to choose
    inline void SetSampleSize(const int s) { m_sampleSize = s; }

    ///Get the priors as a std::array<double, 9>
    inline const std::array<double, 9> GetPriors() const { return m_priors; }
    ///Get the priors by filling a double[9] C array
    inline void GetPriors(double (&p)[9]) const { std::copy(m_priors.begin(), m_priors.end(), std::begin(p)); }
    ///Set the priors from a std::array<double, 9>
    inline void SetPriors(const std::array<double, 9> &p) { m_priors = p; }
    ///Set the priors from a C array
    inline void SetPriors(const double (&p)[9]) { std::copy(std::begin(p), std::end(p), m_priors.begin()); }

protected:
    ///The vectors q1 and q2 are computed from the stain priors and the q adjustment factor
    void ComputeQVectorsFromPriors(cv::InputArray stainPriors, cv::OutputArray qVectors, double qVectorMixRatio);

private:
    double m_avgODThreshold;
    double m_percentileThreshold;
    double m_qVectorMixRatio;

    ///Stain vectors to use as a starting point for the fit
    std::array<double, 9> m_priors;

    ///The number of pixels that should be used to calculate the stain vectors
    int m_sampleSize;
};

} // namespace image
} // namespace sedeen
#endif
