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

//OpenCV include
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>

namespace sedeen {
namespace image {

class MacenkoHistogram {

public:
    MacenkoHistogram();
    ~MacenkoHistogram();

    ///Given a set of 2D vectors (rows), find angle (w/ atan2), histogram, find vectors at hi/lo percentile thresholds
    bool PercentileThresholdVectors(cv::InputArray projectedPoints, cv::OutputArray percentileThreshPoints, 
        const double percentileThresholdValue);
    ///Two-parameter overload of PercentileThresholdVectors, uses the member variable value for the percentileThresholdValue
    bool PercentileThresholdVectors(cv::InputArray projectedPoints, cv::OutputArray percentileThreshPoints);

    ///Convert a set of 2D vectors to float angles between -pi and pi using the arctan2 function
    void VectorsToAngles(cv::InputArray inputVectors, cv::OutputArray outputAngles);
    ///Convert a set of angles to 2D vectors (CV input)
    void AnglesToVectors(cv::InputArray inputAngles, cv::OutputArray outputVectors);
    ///Convert a set of angles to 2D vectors (std::array input)
    void AnglesToVectors(const std::array<double,2> &inputAngles, cv::OutputArray outputVectors);

    ///Create a histogram of angle values, assumed to be between -pi and pi, find angles at %ile thresholds
    const std::array<double, 2> FindPercentileThresholdValues(cv::InputArray vals);

    ///Set the percentileThreshold member variable (force to be between 0 and 50%)
    inline void SetPercentileThreshold(const double &p) { 
        double _p = p < 0.0 ? 0.0 : p; //can't be less than 0
        _p = _p > 100.0 ? 100.0 : _p; //can't be more than 100
        m_percentileThreshold = (_p <= 50.0) ? _p : (100.0 - _p); //100-p if over 50
    }
    ///Get the percentileThreshold member variable
    inline const double GetPercentileThreshold() const { return m_percentileThreshold; }

    ///Set/Get the number of histogram bins
    inline void SetNumHistogramBins(const int n) { m_numHistogramBins = n; }
    ///Set/Get the number of histogram bins
    inline const int GetNumHistogramBins() const { return m_numHistogramBins; }

    ///Set/Get the histogram range
    inline void SetHistogramRange(const std::array<float, 2> &r) { m_histRange = r; }
    ///Set/Get the histogram range
    inline const std::array<float, 2> GetHistogramRange() const { return m_histRange; }

private:
    double m_percentileThreshold;
    int m_numHistogramBins;
    std::array<float, 2> m_histRange;
};

} // namespace image
} // namespace sedeen
#endif
