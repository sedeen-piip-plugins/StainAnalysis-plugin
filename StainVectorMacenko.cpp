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


//I'm moving this, but keeping the code snippets for later use
    //Start the steady clock
    //auto startTime = std::chrono::steady_clock::now();

    //How long did it take to get the basis vectors?
    //auto afterBasisVectorsTime = std::chrono::steady_clock::now();

    //How long did getting the basis vectors take?
    //long long gettingBasisVectorsDuration = std::chrono::duration_cast<std::chrono::microseconds>(afterBasisVectorsTime - startTime).count();



#include "StainVectorMacenko.h"

//#include <chrono>
//#include <random>
#include <fstream>
#include <sstream>

//For now, include ODConversion here, but try to do the OD conversion and thresholding
//in a kernel, and use a factory to apply it before passing the factory to this class
#include "ODConversion.h"
#include "StainVectorMath.h"
#include "MacenkoHistogram.h"
#include "BasisTransform.h"

namespace sedeen {
namespace image {

StainVectorMacenko::StainVectorMacenko(std::shared_ptr<tile::Factory> source)
    : StainVectorBase(source),
    m_sampleSize(0), //Must set to greater than 0 to ComputeStainVectors
    m_avgODThreshold(0.15), //assign default value
    m_percentileThreshold(1.0) //assign default value
{}//end constructor

StainVectorMacenko::~StainVectorMacenko(void) {
}//end destructor

void StainVectorMacenko::ComputeStainVectors(double outputVectors[9]) {
    if (this->GetSourceFactory() == nullptr) { return; }
    //Using this overload of the method requires setting sample size in advance
    int sampleSize = this->GetSampleSize();
    if (sampleSize <= 0) { return; }
    double ODthreshold = this->GetODThreshold();
    double percentileThreshold = this->GetPercentileThreshold();
    if (percentileThreshold <= 0.0) { return; }

    //Sample a set of pixel values from the source
    cv::Mat samplePixels;
    auto theSampler = this->GetRandomWSISampler();
    if (theSampler == nullptr) { return; }
    bool samplingSuccess = theSampler->ChooseRandomPixels(samplePixels, sampleSize, ODthreshold);
    if (!samplingSuccess) { return; }

    //Temp file output
    std::fstream tempOut;
    tempOut.open("D:\\mschumaker\\projects\\Sedeen\\testData\\output\\tempout-ComputeStainVectors.txt", std::fstream::out);


    //Create a class to perform the basis transformation of the sample pixels
    std::unique_ptr<BasisTransform> theBasisTransform = std::make_unique<BasisTransform>();
    //Both the input and output data points should be the matrix rows (columns are pixel elements)
    cv::Mat projectedPoints;
    theBasisTransform->PCAPointTransform(samplePixels, projectedPoints);

    //Now, the only reason to see the basis vectors in this class is out of curiosity
    //cv::Mat basisVectors;
    //bool getBasisSuccess = theBasisTransform->GetBasisVectors(basisVectors);

    std::stringstream scov;
    //scov << "The basis vectors are: " << basisVectors << std::endl;
    //scov << "the projectedPoints are: " << std::endl;
    //scov << projectedPoints << std::endl;




    //Create a class to histogram the results and find 2D vectors corresponding to percentile thresholds
    std::unique_ptr<MacenkoHistogram> theHistogram = std::make_unique<MacenkoHistogram>();
    cv::Mat percentileThreshVectors;
    theHistogram->PercentileThresholdVectors(projectedPoints, percentileThreshVectors, this->GetPercentileThreshold());

    scov << "The percentile threshold vectors: " << std::endl;
    scov << percentileThreshVectors << std::endl;


    //Back-project to get un-normalized stain vectors. DO NOT translate to the mean after backprojection.
    cv::Mat backProjectedVectors;
    theBasisTransform->backProjectPoints(percentileThreshVectors, backProjectedVectors, false); //useMean=false

    scov << "The back projected percentileThreshVectors: " << std::endl;
    scov << backProjectedVectors << std::endl;


    //Back-project to get un-normalized stain vectors
    //cv::Mat stainVectorOutput;


    //Write to a temp array first
    double tempStainVecOutput[9] = {0.0};



    //MAKE THIS GOOD

    StainCVMatToCArray(backProjectedVectors, tempStainVecOutput);
    
    
    
    
    //Error checks?
    for (int i = 0; i < 9; i++) {
        outputVectors[i] = tempStainVecOutput[i];
    }

    tempOut << scov.str() << std::endl;
    tempOut.close();



    //Test converting it the other way
    cv::Mat convertBack;
    StainCArrayToCVMat(tempStainVecOutput, convertBack);



}//end single-parameter ComputeStainVectors



//This overload does not have a default value for sampleSize, so it requires at least two arguments,
//thus there is a clear difference in arguments between this and the other overload of the method
void StainVectorMacenko::ComputeStainVectors(double outputVectors[9], int sampleSize,
    const double ODthreshold /* = 0.15 */, const double percentileThreshold /* = 1.0 */) {
    if (this->GetSourceFactory() == nullptr) { return; }
    //Set member variables with the argument values
    this->SetSampleSize(sampleSize);
    this->SetODThreshold(ODthreshold);
    this->SetPercentileThreshold(percentileThreshold);
    //Call the single-parameter version of this method, which uses the member variables
    this->ComputeStainVectors(outputVectors);
}//end multi-parameter ComputeStainVectors



void StainVectorMacenko::StainCVMatToCArray(cv::InputArray inputData, double outputVectors[9]) {
    if (inputData.empty()) { return; }
    cv::Mat stainVectorOutput, _stainVectorMat(inputData.getMat());
    _stainVectorMat.convertTo(stainVectorOutput, cv::DataType<double>::type);
    //Reshape the matrix and get the data
    int numElements = static_cast<int>(stainVectorOutput.total());


    //pick up here!!!



    //Fill a temporary matrix with the stain vector elements before normalization
    double nonUnitary[9] = { 0.0 };




    //Normalize the stain vectors
    //StainVectorMath::Make3x3MatrixUnitary(nonUnitary, outputVectors);




    //scov << "Well, here are what it calculates as the stain vectors:" << std::endl;
    //scov << stainVectorOutput << std::endl;

    //double rowSums[2] = { 0.0 };
    //rowSums[0] = stainVectorOutput.at<double>(0, 0)
    //    + stainVectorOutput.at<double>(0, 1)
    //    + stainVectorOutput.at<double>(0, 2);
    //rowSums[1] = stainVectorOutput.at<double>(1, 0)
    //    + stainVectorOutput.at<double>(1, 1)
    //    + stainVectorOutput.at<double>(1, 2);

    //if (rowSums[0] < 0.0) {
    //    stainVectorOutput.at<double>(0, 0) *= -1.0;
    //    stainVectorOutput.at<double>(0, 1) *= -1.0;
    //    stainVectorOutput.at<double>(0, 2) *= -1.0;
    //}
    //if (rowSums[1] < 0.0) {
    //    stainVectorOutput.at<double>(1, 0) *= -1.0;
    //    stainVectorOutput.at<double>(1, 1) *= -1.0;
    //    stainVectorOutput.at<double>(1, 2) *= -1.0;
    //}


    //scov << "Here they are positive summed:" << std::endl;
    //scov << stainVectorOutput << std::endl;


    nonUnitary[0] = stainVectorOutput.at<double>(0, 0);
    nonUnitary[1] = stainVectorOutput.at<double>(0, 1);
    nonUnitary[2] = stainVectorOutput.at<double>(0, 2);
    nonUnitary[3] = stainVectorOutput.at<double>(1, 0);
    nonUnitary[4] = stainVectorOutput.at<double>(1, 1);
    nonUnitary[5] = stainVectorOutput.at<double>(1, 2);
    nonUnitary[6] = 0.0;
    nonUnitary[7] = 0.0;
    nonUnitary[8] = 0.0;

    StainVectorMath::Make3x3MatrixUnitary(nonUnitary, outputVectors);
    //scov << "Here they are unitary: " << std::endl;

    //scov << outputVectors[0] << ", " << outputVectors[1] << ", " << outputVectors[2] << std::endl;
    //scov << outputVectors[3] << ", " << outputVectors[4] << ", " << outputVectors[5] << std::endl;
    //scov << outputVectors[6] << ", " << outputVectors[7] << ", " << outputVectors[8] << std::endl;






 /*   scov << "Well, here are what it calculates as the stain vectors:" << std::endl;
    scov << stainVectorOutput << std::endl;

    double rowSums[2] = { 0.0 };
    rowSums[0] = stainVectorOutput.at<double>(0, 0)
        + stainVectorOutput.at<double>(0, 1)
        + stainVectorOutput.at<double>(0, 2);
    rowSums[1] = stainVectorOutput.at<double>(1, 0)
        + stainVectorOutput.at<double>(1, 1)
        + stainVectorOutput.at<double>(1, 2);

    if (rowSums[0] < 0.0) {
        stainVectorOutput.at<double>(0, 0) *= -1.0;
        stainVectorOutput.at<double>(0, 1) *= -1.0;
        stainVectorOutput.at<double>(0, 2) *= -1.0;
    }
    if (rowSums[1] < 0.0) {
        stainVectorOutput.at<double>(1, 0) *= -1.0;
        stainVectorOutput.at<double>(1, 1) *= -1.0;
        stainVectorOutput.at<double>(1, 2) *= -1.0;
    }


    scov << "Here they are positive summed:" << std::endl;
    scov << stainVectorOutput << std::endl;

    double nonUnitary[9] = { 0.0 };
    nonUnitary[0] = stainVectorOutput.at<double>(0, 0);
    nonUnitary[1] = stainVectorOutput.at<double>(0, 1);
    nonUnitary[2] = stainVectorOutput.at<double>(0, 2);
    nonUnitary[3] = stainVectorOutput.at<double>(1, 0);
    nonUnitary[4] = stainVectorOutput.at<double>(1, 1);
    nonUnitary[5] = stainVectorOutput.at<double>(1, 2);
    nonUnitary[6] = 0.0;
    nonUnitary[7] = 0.0;
    nonUnitary[8] = 0.0;

    scov << "Here they are unitary: " << std::endl;

    scov << outputVectors[0] << ", " << outputVectors[1] << ", " << outputVectors[2] << std::endl;
    scov << outputVectors[3] << ", " << outputVectors[4] << ", " << outputVectors[5] << std::endl;
    scov << outputVectors[6] << ", " << outputVectors[7] << ", " << outputVectors[8] << std::endl;
*/




}//end StainCVMatToCArray



void StainVectorMacenko::StainCArrayToCVMat(double inputVectors[9], cv::OutputArray outputData) {
    //Create a copy of the input data
    double inputCopy[9] = { 0.0 };
    for (int i = 0; i < 9; i++) {
        inputCopy[i] = inputVectors[i];
    }
    //Create a cv::Mat of type double and point to inputCopy (reference, not deep copy)
    cv::Mat inputMatFlat(1, 9, cv::DataType<double>::type, inputCopy);
    //Reshape the matrix (does not reallocate, still points to inputCopy)
    cv::Mat inputMatSquare = inputMatFlat.reshape(0, 3);
    outputData.assign(inputMatSquare);
}//end StainCArrayToCVMat

} // namespace image
} // namespace sedeen
