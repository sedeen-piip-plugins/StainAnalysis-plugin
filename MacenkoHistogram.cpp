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

#include "MacenkoHistogram.h"

#include <cmath>
#include <fstream>
#include <sstream>

namespace sedeen {
namespace image {

MacenkoHistogram::MacenkoHistogram() : 
    m_percentileThreshold(1.0),
    m_numHistogramBins(1024),
    m_histRange({ -static_cast<float>(CV_PI),static_cast<float>(CV_PI) })
{
}//end constructor

MacenkoHistogram::~MacenkoHistogram(void) {
}//end destructor


bool MacenkoHistogram::AssignClusters(cv::InputArray projectedPoints, cv::InputOutputArray clusterAssignments, 
    cv::InputArray qPriors) {



    //return true on success
    return true;
}//end AssignClusters





bool MacenkoHistogram::PercentileThresholdVectors(cv::InputArray projectedPoints,
    cv::OutputArray percentileThreshPoints, const double percentileThresholdValue) {
    //SetPercentileThreshold forces the value to be between 0 and 50
    this->SetPercentileThreshold(percentileThresholdValue);
    //Call the 2-parameter overload of this method,
    //after the member variable value has been set
    return this->PercentileThresholdVectors(projectedPoints, percentileThreshPoints);
}//end PercentileThresholdVectors

bool MacenkoHistogram::PercentileThresholdVectors(cv::InputArray projectedPoints, 
    cv::OutputArray percentileThreshPoints) {
    //Check the value of the member variable, return false if it is out of range
    double threshVal = this->GetPercentileThreshold();
    if ((threshVal <= 0.0) || (threshVal >= 100.0)) {
        return false;
    }

    //Get the angular coordinates of the 2D points
    cv::Mat angleVals;
    this->VectorsToAngles(projectedPoints, angleVals);
    //Check if angleVals is empty, return failure if so
    if (angleVals.empty()) { return false; }

    //Histogram the angles, get angles at percentile values
    std::array<double,2> percentileAngles = FindPercentileThresholdValues(angleVals);
    if (percentileAngles.empty()) { return false; }

    cv::Mat angToVecOutput;
    this->AnglesToVectors(percentileAngles, angToVecOutput);
    if (angToVecOutput.empty()) { return false; }

    //Return true on success
    percentileThreshPoints.assign(angToVecOutput);
    return true;
}//end PercentileThresholdVectors

void MacenkoHistogram::VectorsToAngles(cv::InputArray inputVectors, cv::OutputArray outputAngles) {
    //Check if inputVectors not empty, and there are at least two columns
    if (inputVectors.empty()) { return; }
    if (inputVectors.cols() < 2) { return; }

    //Access the input data as a matrix
    cv::Mat inputMat = inputVectors.getMat();
    cv::Mat inputFloatMat;
    inputMat.convertTo(inputFloatMat, cv::DataType<float>::type);

    //Define the output matrix
    cv::Mat angleVals(inputVectors.rows(), 1, cv::DataType<float>::type);
    //Calculate output matrix values
    for (auto p = angleVals.begin<float>(); p != angleVals.end<float>(); p++) {
        int row = static_cast<int>(p.lpos());
        //atan2 is undefined when x and y are both 0
        bool angleUndef = (std::abs(inputFloatMat.at<float>(row, 0)) < 1e-4) 
            && (std::abs(inputFloatMat.at<float>(row, 1)) < 1e-4); //Choice of zero roundoff threshold is arbitrary
        //Use maximum float value as the undefined value
        float undefined = std::numeric_limits<float>::max();
        //Calculate the arctan2, unless undefined
        *p = angleUndef ? undefined
            : std::atan2(inputFloatMat.at<float>(row, 1), inputFloatMat.at<float>(row, 0));
    }
    outputAngles.assign(angleVals);
}//end VectorsToAngles

void MacenkoHistogram::AnglesToVectors(cv::InputArray inputAngles, cv::OutputArray outputVectors) {
    //Check if inputAngles is empty
    if (inputAngles.empty()) { return; }
    //Access the input angle data as a matrix
    cv::Mat inputMat = inputAngles.getMat();
    cv::Mat inputDoubleMat;
    inputMat.convertTo(inputDoubleMat, cv::DataType<double>::type);
    std::array<double, 2> angleArray;
    if (inputDoubleMat.rows == 1 && inputDoubleMat.cols >= 2) {
        angleArray[0] = inputDoubleMat.at<double>(0, 0);
        angleArray[1] = inputDoubleMat.at<double>(0, 1);
    }
    else if (inputMat.rows >= 2 && inputMat.cols == 1) {
        angleArray[0] = inputDoubleMat.at<double>(0, 0);
        angleArray[1] = inputDoubleMat.at<double>(1, 0);
    }
    else {
        return;
    }
    //Call the other overload of this method
    this->AnglesToVectors(angleArray, outputVectors);
}//end AnglesToVectors

void MacenkoHistogram::AnglesToVectors(const std::array<double,2> &inputAngles, cv::OutputArray outputVectors) {
    //Check if inputAngles is empty
    if (inputAngles.empty()) { return; }
    //Check if inputAngles is 0,0
    if (inputAngles[0] == 0.0 && inputAngles[1] == 0.0) { return; }

    //Convert to 2D Cartesian
    cv::Mat cartesianVectors(2, 2, cv::DataType<double>::type);
    cartesianVectors.at<double>(0, 0) = std::cos(inputAngles[0]);
    cartesianVectors.at<double>(0, 1) = std::sin(inputAngles[0]);
    cartesianVectors.at<double>(1, 0) = std::cos(inputAngles[1]);
    cartesianVectors.at<double>(1, 1) = std::sin(inputAngles[1]);
    outputVectors.assign(cartesianVectors);
}//end AnglesToVectors

const std::array<double, 2> MacenkoHistogram::FindPercentileThresholdValues(cv::InputArray vals) {
    //Return threshold angle values as a 2-element array
    //The intent is to make the output format more strict than cv::Mat would be
    std::array<double, 2> errorValues = { 0.0,0.0 };
    std::array<double, 2> outputValues = { 0.0,0.0 };
    if (vals.empty()) { return errorValues; }

    //Temp file output
    std::fstream tempOut;
    tempOut.open("D:\\mschumaker\\projects\\Sedeen\\testData\\output\\tempout-theHistogram.txt", std::fstream::out);
    std::stringstream sh;


    //cv::calcHist only operates on floats, so convert data type
    cv::Mat _vals, floatVals;
    _vals = vals.getMat(); // .clone();

    sh << "The values before conversion: " << std::endl;
    sh << vals.getMat() << std::endl;

    _vals.convertTo(floatVals, cv::DataType<float>::type);

    //Make sure that the histogram range is valid
    std::array<float, 2> rangeArray = this->GetHistogramRange();
    if (rangeArray[1] <= rangeArray[0]) { return errorValues; }
    //Check the current number set for histogram bins
    if (this->GetNumHistogramBins() <= 0) { return errorValues; }

    //Fill the histogram
    cv::Mat theHist;
    int channels[1] = { 0 };
    int histSize[1] = { this->GetNumHistogramBins() };
    float range[] = { rangeArray[0], rangeArray[1] };
    const float* histRange[] = { range };
    bool uniform = true;
    bool accumulate = false;
    cv::calcHist(&floatVals, 1, channels, cv::Mat(), theHist, 1, histSize, histRange, uniform, accumulate);
    //The resulting type of theHist elements is float



    sh << "Here are all of the angles being added to the histogram: " << std::endl;
    sh << floatVals << std::endl;

    sh << "And here's the histogram: " << std::endl;
    sh << theHist << std::endl;


    //Count how many elements were added to the histogram
    cv::Scalar scalarTotal = cv::sum(theHist);
    double histoCountTotal = scalarTotal[0];

    sh << "The value of histoCountTotal is : " << histoCountTotal << std::endl;

    //Calculate the slope and intercept to convert from bin to value
    //value = intercept + slope*bin
    double intercept = static_cast<double>(range[0]);
    double slope = static_cast<double>(range[1] - range[0]) / static_cast<double>(histSize[0]);
    double percentileThreshold = this->GetPercentileThreshold();

    //Find the bins containing the lower and upper percentiles
    double lowerFraction = (percentileThreshold) / 100.0;
    double upperFraction = (100.0 - percentileThreshold) / 100.0;
    bool lowerPassed(false), upperPassed(false);
    double lowerBin(-1.0), upperBin(-1.0);
    double cumulativeSum(0.0);

    sh << "Looking for the lower and upper bins! Time to search! " << std::endl;

    for (int bin = 0; bin < (theHist.size)[0]; bin++) {
        double prevFraction = cumulativeSum / histoCountTotal;
        cumulativeSum += static_cast<double>(theHist.at<float>(bin, 0));
        double currentFraction = cumulativeSum / histoCountTotal;

        sh << "bin =" << bin << ", prevFraction= " << prevFraction << ", cumulativeSum= " << cumulativeSum 
            << ", currentFraction= " << currentFraction << std::endl;

        if (!lowerPassed && (currentFraction >= lowerFraction)) {
            lowerPassed = true;
            //linearly interpolate a bin position

            sh << "This is the bin for which I think the lowerBin was passed! " << bin << std::endl;
            lowerBin = static_cast<double>(bin - 1) + (lowerFraction - prevFraction) / (currentFraction - prevFraction);
            sh << "And here's the value of lowerBin I calculate: " << lowerBin << std::endl;
        }
        if (!upperPassed && (currentFraction >= upperFraction)) {
            upperPassed = true;
            //linearly interpolate a bin position
            sh << "This is the bin for which I think the upperBin was passed! " << bin << std::endl;
            upperBin = static_cast<double>(bin - 1) + (upperFraction - prevFraction) / (currentFraction - prevFraction);
            sh << "And here's the value of upperBin I calculate: " << upperBin << std::endl;
        }
        if (lowerPassed && upperPassed) { break; }
    }
    //Convert from bins to values
    outputValues[0] = intercept + slope * lowerBin;
    outputValues[1] = intercept + slope * upperBin;

    sh << "This is what I think intercept and slope are: " << intercept << " and " << slope << std::endl;
    sh << "and here are the output values I'm returning: " << outputValues[0] << " and " << outputValues[1] << std::endl;

    tempOut << sh.str() << std::endl;
    tempOut.close();


    return outputValues;
}//end FindPercentileThresholdValues

} // namespace image
} // namespace sedeen
