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

#ifndef STAINANALYSIS_NIETHAMMERHISTOGRAM_H
#define STAINANALYSIS_NIETHAMMERHISTOGRAM_H

#include <array>

#include "AngleHistogram.h"

//OpenCV include
#include <opencv2/core/core.hpp>

namespace sedeen {
namespace image {

class NiethammerHistogram : public AngleHistogram {

public:
    ///Constructor with two parameters, with defaults specified
    NiethammerHistogram(double alpha = 0.15, int nbins = 128);
    ///Destructor
    virtual ~NiethammerHistogram();


    //There will be refactoring to separate histogramming to get percentile values from
    //histogramming to find an Otsu threshold.

    ///Computes angles, finds Otsu threshold, reassigns points belonging to clusters (below or above threshold).
    bool AssignClusters(cv::InputArray projectedPoints, cv::InputOutputArray clusterAssignments, cv::InputArray stainPriors);
    
public:
    ///Set/Get the mixing parameter between the raw stain priors (basis vectors) to get the mixed vectors q1 and q2
    inline void SetAlphaMixRatio(const double &alpha) { m_alphaMixRatio = alpha; }
    ///Set/Get the mixing parameter between the raw stain priors (basis vectors) to get the mixed vectors q1 and q2
    inline const double GetAlphaMixRatio() const { return m_alphaMixRatio; }



private:
    ///Parameter that determines the amount of mixing between the raw stain priors (basis vectors)
    double m_alphaMixRatio;
    ///The raw stain priors, as used externally
    cv::Mat m_stainPriors;
    ///The adjusted priors, calculated by mixing the raw stain priors using the value of m_alphaMixRatio
    cv::Mat m_qPriors;

};

} // namespace image
} // namespace sedeen
#endif

