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

#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>

namespace sedeen {
namespace image {

MacenkoHistogram::MacenkoHistogram() : m_numTestingPixels(10),
    m_rgen((std::random_device())()) //Initialize random number generation
{
}//end constructor

MacenkoHistogram::~MacenkoHistogram(void) {
}//end destructor

void MacenkoHistogram::OptimizeBasisVectorSigns(cv::InputArray sourcePixels, /*assume sourcePixels to be row vectors */
    cv::InputArray inputVectors, cv::OutputArray outputVectors,
    VectorDirection basisVecDir /* = VectorDirection::COLUMNVECTORS */) {
    //Check a small number of source pixels to try to get projected points with all positive elements
    //The number of pixels to check is arbitrary, default value 10
    cv::Mat subsampleofPixels;
    int numTestPixels = this->GetNumTestingPixels();
    if (numTestPixels < 1) {
        //Deep copy the inputVectors to the outputVectors
        outputVectors.assign(inputVectors.getMat().clone());
        return;
    }
    else {
        CreatePixelSubsample(sourcePixels, subsampleofPixels, numTestPixels);
    }

    //Arrange the basis vectors as the columns, if they're not already
    cv::Mat tempVectors;
    if (basisVecDir == VectorDirection::COLUMNVECTORS) {
        //no change
        tempVectors = inputVectors.getMat().clone();
    }
    else if (basisVecDir == VectorDirection::ROWVECTORS) {
        //transpose the inputVectors matrix
        cv::transpose(inputVectors, tempVectors);
    }
    else {
        //invalid, but use the inputVectors as given and proceed anyway
        tempVectors = inputVectors.getMat().clone();
    }


    //Temp file output
    std::fstream tempOut;
    tempOut.open("D:\\mschumaker\\projects\\Sedeen\\testData\\output\\tempout-samplepix.txt", std::fstream::out);
    std::stringstream ss;


    



    //Create the list of +/- options for the basis vectors (0 is +, 1 is -)
    int numCombinations = static_cast<int>(std::pow(2, tempVectors.cols));
    cv::Mat testMultFactors(numCombinations, tempVectors.cols, cv::DataType<int>::type);
    for (int row = 0; row < numCombinations; row++) {
        for (int col = 0; col < tempVectors.cols; col++) {
            int newVal = (row >> col) & 1; //bit shift and mask
            testMultFactors.at<int>(row, col) = newVal;
        }
    }

    ss << testMultFactors << std::endl;



    //Now do something with the subsampleOfPixels




    tempOut << ss.str() << std::endl;
    tempOut.close();


}//end OptimizeBasisVectorSigns



void MacenkoHistogram::CreatePixelSubsample(cv::InputArray sourcePixels, cv::OutputArray subsample, int numberOfPixels) {
    cv::Mat tempSubsampleMat;
    if (numberOfPixels < 1) {
        return;
    }
    else if (numberOfPixels <= sourcePixels.rows()) {
        //Use all of the source pixels
        subsample.assign(sourcePixels.getMat().clone());
    }
    else {
        //Create with no rows, push new rows later
        tempSubsampleMat = cv::Mat(0, sourcePixels.cols(), cv::DataType<double>::type);
    }

    //Initialize a uniform distribution to choose random pixels
    std::uniform_int_distribution<int> randSourcePixelIndex(0, sourcePixels.rows() - 1);
    //Create a list of pixel (row) indices for the randomized subset, without duplication
    std::vector<int> pixelList;
    for (int px = 0; px < numberOfPixels; px++) {
        int countLimit = 2 * numberOfPixels; //loop count limit (just in case)
        bool newIndexFound = false;
        int attemptNumber = 0;
        while (!newIndexFound && (attemptNumber < countLimit)) {
            int newIndex = randSourcePixelIndex(m_rgen);
            auto it = std::find(pixelList.begin(), pixelList.end(), newIndex);
            if (it == pixelList.end()) {
                newIndexFound = true;
                pixelList.push_back(newIndex);
            }
            else {
                newIndexFound = false;
                attemptNumber++;
            }
        }
    }
    //Get pixels at the row indices in pixelList, push to tempSubsampleMat
    for (auto pxit = pixelList.begin(); pxit != pixelList.end(); ++pxit) {
        tempSubsampleMat.push_back(sourcePixels.getMat().row(*pxit));
    }
    //assign to the OutputArray
    subsample.assign(tempSubsampleMat);
}//end CreatePixelSubsample



} // namespace image
} // namespace sedeen
