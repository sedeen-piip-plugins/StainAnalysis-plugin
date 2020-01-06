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

#include "MacenkoHistogram.h"

#include <cmath>
#include <fstream>
#include <sstream>

namespace sedeen {
namespace image {

MacenkoHistogram::MacenkoHistogram(double pthresh /* = 1.0 */, int nbins /*= 1024 */) : 
    AngleHistogram(nbins /*, default range */), m_percentileThreshold(pthresh) {
}//end constructor

MacenkoHistogram::~MacenkoHistogram(void) {
}//end destructor

bool MacenkoHistogram::PercentileThresholdVectors(cv::InputArray projectedPoints,
    cv::OutputArray percentileThreshPoints, const double &percentileThresholdValue) {
    //SetPercentileThreshold forces the value to be between 0 and 50
    this->SetPercentileThreshold(percentileThresholdValue);
    //Call the 2-parameter overload of this method,
    //after the member variable value has been set
    return this->PercentileThresholdVectors(projectedPoints, percentileThreshPoints);
}//end PercentileThresholdVectors

bool MacenkoHistogram::PercentileThresholdVectors(cv::InputArray projectedPoints, 
    cv::OutputArray percentileThreshPoints) {
    //Check the value of the member variable, return false if it is out of range
    float threshVal = static_cast<float>(this->GetPercentileThreshold());
    if ((threshVal <= 0.0) || (threshVal >= 100.0)) {
        return false;
    }

    //Get the angular coordinates of the 2D points
    cv::Mat angleVals;
    this->VectorsToAngles(projectedPoints, angleVals);
    //Check if angleVals is empty, return failure if so
    if (angleVals.empty()) { return false; }

    //Histogram the angles
    cv::Mat theAngleHist;
    FillHistogram(angleVals, theAngleHist);

    std::array<float,2> percentileAngles = FindPercentileThresholdValues(theAngleHist);
    if (percentileAngles.empty()) { return false; }

    cv::Mat angToVecOutput;
    this->AnglesToVectors(percentileAngles, angToVecOutput);
    if (angToVecOutput.empty()) { return false; }

    //Return true on success
    percentileThreshPoints.assign(angToVecOutput);
    return true;
}//end PercentileThresholdVectors



const std::array<float, 2> MacenkoHistogram::FindPercentileThresholdValues(cv::InputArray _theHist) {
    //Return percentile threshold values in the histogram as a 2-element array
    //The intent is to make the output format more strict than cv::Mat would be
    std::array<float, 2> errorValues = { 0.0,0.0 };
    std::array<float, 2> outputValues = { 0.0,0.0 };
    if (_theHist.empty()) { return errorValues; }
    //Get the histogram as a float matrix
    cv::Mat _hist = _theHist.getMat();
    cv::Mat theHistMat;
    _hist.convertTo(theHistMat, cv::DataType<float>::type);

    //Temp file output
    std::fstream tempOut;
    tempOut.open("D:\\mschumaker\\projects\\Sedeen\\testData\\output\\tempout-percentileThreshValues.txt", std::fstream::out);
    std::stringstream sh;

    sh << "Here's the histogram: " << std::endl;
    sh << theHistMat << std::endl;
    
    //Count how many elements were added to the histogram
    cv::Scalar scalarTotal = cv::sum(theHistMat);
    float histoCountTotal = scalarTotal[0];

    sh << "The value of histoCountTotal is : " << histoCountTotal << std::endl;

    //Find the bins containing the lower and upper percentiles
    float percentileThreshold = static_cast<float>(this->GetPercentileThreshold());
    float lowerFraction = (percentileThreshold) / 100.0;
    float upperFraction = (100.0 - percentileThreshold) / 100.0;
    bool lowerPassed(false), upperPassed(false);
    float lowerBin(-1.0), upperBin(-1.0);
    float cumulativeSum(0.0);

    sh << "Looking for the lower and upper bins! Time to search! " << std::endl;

    for (int bin = 0; bin < (theHistMat.size)[0]; bin++) {
        float prevFraction = cumulativeSum / histoCountTotal;
        cumulativeSum += static_cast<float>(theHistMat.at<float>(bin, 0));
        float currentFraction = cumulativeSum / histoCountTotal;

        sh << "bin =" << bin << ", prevFraction= " << prevFraction << ", cumulativeSum= " << cumulativeSum 
            << ", currentFraction= " << currentFraction << std::endl;

        if (!lowerPassed && (currentFraction >= lowerFraction)) {
            lowerPassed = true;
            //linearly interpolate a bin position

            sh << "This is the bin for which I think the lowerBin was passed! " << bin << std::endl;
            lowerBin = static_cast<float>(bin - 1) + (lowerFraction - prevFraction) / (currentFraction - prevFraction);
            sh << "And here's the value of lowerBin I calculate: " << lowerBin << std::endl;
        }
        if (!upperPassed && (currentFraction >= upperFraction)) {
            upperPassed = true;
            //linearly interpolate a bin position
            sh << "This is the bin for which I think the upperBin was passed! " << bin << std::endl;
            upperBin = static_cast<float>(bin - 1) + (upperFraction - prevFraction) / (currentFraction - prevFraction);
            sh << "And here's the value of upperBin I calculate: " << upperBin << std::endl;
        }
        if (lowerPassed && upperPassed) { break; }
    }
    //Convert from bins to values
    outputValues[0] = HistogramBinToAngle(lowerBin);
    outputValues[1] = HistogramBinToAngle(upperBin);

    //sh << "This is what I think intercept and slope are: " << intercept << " and " << slope << std::endl;
    sh << "and here are the output values I'm returning: " << outputValues[0] << " and " << outputValues[1] << std::endl;

    tempOut << sh.str() << std::endl;
    tempOut.close();


    return outputValues;
}//end FindPercentileThresholdValues

} // namespace image
} // namespace sedeen
