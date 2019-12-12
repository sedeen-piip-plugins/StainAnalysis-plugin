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

#include "AngleHistogram.h"

#include <cmath>
#include <fstream>
#include <sstream>

namespace sedeen {
namespace image {

AngleHistogram::AngleHistogram(int nbins /*= 128 */, std::array<float, 2> range /* -CV_PI to CV_PI */) :
    m_numHistogramBins(nbins), m_histRange(range) {
}//end constructor

AngleHistogram::~AngleHistogram(void) {
}//end destructor

void AngleHistogram::FillHistogram(cv::InputArray inVals, cv::OutputArray outHist) {
    //Get the member variables for number of bins and range, pass to 4-parameter overload of this method
    if (inVals.empty()) { return; }
    //Make sure that the histogram range is valid
    std::array<float, 2> rangeArray = this->GetHistogramRange();
    if (rangeArray[1] <= rangeArray[0]) { return; }
    //Check the current number set for histogram bins
    int nbins = this->GetNumHistogramBins();
    if (nbins <= 0) { return; }
    //Invoke the 4-parameter overload
    FillHistogram(inVals, outHist, nbins, rangeArray);
}//end FillHistogram (public 2-parameter)


void AngleHistogram::FillHistogram(cv::InputArray inVals, cv::OutputArray outHist,
    int nbins, std::array<float, 2> rangeArray) {
    if (inVals.empty()) { return; }
    if (rangeArray[1] <= rangeArray[0]) { return; }
    if (nbins <= 0) { return; }

    //Temp file output
    std::fstream tempOut;
    tempOut.open("D:\\mschumaker\\projects\\Sedeen\\testData\\output\\tempout-theHistogram.txt", std::fstream::out);
    std::stringstream sh;

    //cv::calcHist only operates on floats, so convert data type
    cv::Mat _vals, floatVals;
    _vals = inVals.getMat();

    sh << "The values before conversion: " << std::endl;
    sh << inVals.getMat() << std::endl;

    _vals.convertTo(floatVals, cv::DataType<float>::type);

    //Fill the histogram
    cv::Mat theHist;
    int channels[1] = { 0 };
    int histSize[1] = { nbins };
    float range[] = { rangeArray[0], rangeArray[1] };
    const float* histRange[] = { range };
    bool uniform = true;
    bool accumulate = false;
    cv::calcHist(&floatVals, 1, channels, cv::Mat(), theHist, 1, histSize, histRange, uniform, accumulate);
    //The resulting type of theHist elements is float

    sh << "Here's the histogram: " << std::endl;
    sh << theHist << std::endl;

    tempOut << sh.str() << std::endl;
    tempOut.close();

    outHist.assign(theHist);
}//end FillHistogram (protected 4-parameter overload)


void AngleHistogram::VectorsToAngles(cv::InputArray inputVectors, cv::OutputArray outputAngles) {
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
        bool angleUndef = (std::abs(inputFloatMat.at<float>(row, 0)) < 1e-6)
            && (std::abs(inputFloatMat.at<float>(row, 1)) < 1e-6); //Choice of zero roundoff threshold is arbitrary
        //Use maximum float value as the undefined value
        float undefined = std::numeric_limits<float>::max();
        //Calculate the arctan2, unless undefined
        *p = angleUndef ? undefined
            : std::atan2(inputFloatMat.at<float>(row, 1), inputFloatMat.at<float>(row, 0));
    }
    outputAngles.assign(angleVals);
}//end VectorsToAngles

void AngleHistogram::AnglesToVectors(cv::InputArray inputAngles, cv::OutputArray outputVectors) {
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

void AngleHistogram::AnglesToVectors(const std::array<double, 2> &inputAngles, cv::OutputArray outputVectors) {
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

} // namespace image
} // namespace sedeen