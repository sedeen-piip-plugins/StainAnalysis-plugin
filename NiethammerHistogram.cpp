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

#include "NiethammerHistogram.h"

#include <cmath>
#include <fstream>
#include <sstream>

namespace sedeen {
namespace image {

NiethammerHistogram::NiethammerHistogram(double alpha /*= 0.15 */, int nbins /*= 128 */) : AngleHistogram(nbins /*, default range */),
    m_alphaMixRatio(alpha) {
}//end constructor

NiethammerHistogram::~NiethammerHistogram(void) {
}//end destructor


bool NiethammerHistogram::AssignClusters(cv::InputArray projectedPoints, cv::InputOutputArray clusterAssignments,
    cv::InputArray qPriors) {
    //Get the angular coordinates of the 2D points
    cv::Mat angleVals;
    this->VectorsToAngles(projectedPoints, angleVals);
    //Check if angleVals is empty, return failure if so
    if (angleVals.empty()) { return false; }

    //Fill the histogram
    cv::Mat theHist;
    FillHistogram(angleVals, theHist); //Use the member variables for number of bins and histogram range


    //Now what do I do with it?
    //I think I'm right. I test every bin and see which one has the lowest energy function
    //I'm going to need a temporary cluster assignment matrix
    cv::Mat tempClusterAssignments(angleVals.size(), cv::DataType<int>::type);
    //Order: set I_theta to a bin value. 

    //I will have to convert back from histogram bin to angle





    //return true on success
    return true;
}//end AssignClusters




} // namespace image
} // namespace sedeen
