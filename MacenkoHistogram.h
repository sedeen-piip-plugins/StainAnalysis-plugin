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

#ifndef STAINANALYSIS_MACENKOHISTOGRAM_H
#define STAINANALYSIS_MACENKOHISTOGRAM_H

#include <array>

#include "AngleHistogram.h"

//OpenCV include
#include <opencv2/core/core.hpp>

namespace sedeen {
namespace image {

class MacenkoHistogram : public AngleHistogram {

public:
    ///Constructor with two parameters, with defaults specified
    MacenkoHistogram(double pthresh = 1.0, int nbins = 1024);
    ///Destructor
    virtual ~MacenkoHistogram();

    ///Given a set of 2D vectors (rows), find angle (w/ atan2), histogram, find vectors at hi/lo percentile thresholds
    bool PercentileThresholdVectors(cv::InputArray projectedPoints, cv::OutputArray percentileThreshPoints, 
        const double &percentileThresholdValue);

    ///Two-parameter overload of PercentileThresholdVectors, uses the member variable value for the percentileThresholdValue
    bool PercentileThresholdVectors(cv::InputArray projectedPoints, cv::OutputArray percentileThreshPoints);

    ///Given a histogram with range and nbins set in member variables, find values at percentile thresholds
    const std::array<float, 2> FindPercentileThresholdValues(cv::InputArray theHist);

public:
    ///Set the percentileThreshold member variable (force to be between 0 and 50%)
    inline void SetPercentileThreshold(const double &_p) { 
        double p = _p < 0.0 ? 0.0 : _p; //can't be less than 0
        p = p > 100.0 ? 100.0 : p;      //can't be more than 100
        m_percentileThreshold = (p <= 50.0) ? p : (100.0 - p); //100-p if over 50
    }
    ///Get the percentileThreshold member variable
    inline const double GetPercentileThreshold() const { return m_percentileThreshold; }

private:
    double m_percentileThreshold;
};

} // namespace image
} // namespace sedeen
#endif
