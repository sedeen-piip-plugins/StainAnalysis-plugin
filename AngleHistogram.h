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

#ifndef STAINANALYSIS_ANGLEHISTOGRAM_H
#define STAINANALYSIS_ANGLEHISTOGRAM_H

#include <array>

 //OpenCV include
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>

namespace sedeen {
namespace image {

class AngleHistogram {

public:
    ///Constructor with two parameters with defaults specified
    AngleHistogram(int nbins = 128, std::array<float, 2> range = { -static_cast<float>(CV_PI),static_cast<float>(CV_PI) } );
    ///Destructor
    virtual ~AngleHistogram(void);

    ///Populate a histogram from an input array of single-column data, get histogram configuration from member variables
    void FillHistogram(cv::InputArray inVals, cv::OutputArray outHist);

public:
    ///Convert a set of 2D vectors to float angles between -pi and pi using the arctan2 function
    void VectorsToAngles(cv::InputArray inputVectors, cv::OutputArray outputAngles);
    ///Convert a set of angles to 2D vectors (CV input)
    void AnglesToVectors(cv::InputArray inputAngles, cv::OutputArray outputVectors);
    ///Convert a set of angles to 2D vectors (std::array input)
    void AnglesToVectors(const std::array<double, 2> &inputAngles, cv::OutputArray outputVectors);

public:
    ///Set/Get the number of histogram bins
    inline void SetNumHistogramBins(const int &n) { m_numHistogramBins = n; }
    ///Set/Get the number of histogram bins
    inline const int GetNumHistogramBins() const { return m_numHistogramBins; }

    ///Set/Get the histogram range
    inline void SetHistogramRange(const std::array<float, 2> &r) { m_histRange = r; }
    ///Set/Get the histogram range
    inline const std::array<float, 2> GetHistogramRange() const { return m_histRange; }

protected:
    ///Populate a histogram from an input array of single-column data, histogram configuration set by 3rd and 4th arguments
    void FillHistogram(cv::InputArray inVals, cv::OutputArray outHist, int nbins, std::array<float, 2> range);


    //I think I need a way to convert from bin back to angle

    //AngleToHistogramBin

    //HistogramBinToAngle


    //interpolation?




private:
    int m_numHistogramBins;
    std::array<float, 2> m_histRange;
};

} // namespace image
} // namespace sedeen
#endif

