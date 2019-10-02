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

#include "StainVectorMacenko.h"

#include <chrono>
#include <random>
#include <fstream>
#include <sstream>

 //For now, include ODConversion here, but try to do the OD conversion and thresholding
 //in a kernel, and use a factory to apply it before passing the factory to this class
#include "ODConversion.h"
#include "StainVectorMath.h"
#include "MacenkoHistogram.h"

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
    assert(nullptr != this->GetSourceFactory());
    //Using this overload of the method requires setting sample size in advance
    int sampleSize = this->GetSampleSize();
    if (sampleSize <= 0) { return; }
    double ODthreshold = this->GetODThreshold();
    double percentileThreshold = this->GetPercentileThreshold();
    if (percentileThreshold <= 0.0) { return; }

    //Sample a set of pixel values from the source
    cv::Mat samplePixels;
    auto theSampler = this->GetRandomWSISampler();
    assert(nullptr != theSampler);
    bool samplingSuccess = theSampler->ChooseRandomPixels(samplePixels, sampleSize, ODthreshold);
    if (!samplingSuccess) { return; }



    //Temp file output
    std::fstream tempOut;
    tempOut.open("D:\\mschumaker\\projects\\Sedeen\\testData\\output\\tempout.txt", std::fstream::out);
    //std::stringstream sv;

    //sv << samplePixels << std::endl;


    //Let's write these to a file to check that something's happening
    //sv << numPixelsAddedToMatrix-1 << ") Tile " << tl << ", Pixel " << px << ": ";
    //sv << sampledPixelsMatrix.at<double>(numPixelsAddedToMatrix-1, 0) << ", ";
    //sv << sampledPixelsMatrix.at<double>(numPixelsAddedToMatrix-1, 1) << ", ";
    //sv << sampledPixelsMatrix.at<double>(numPixelsAddedToMatrix-1, 2);
    //tempOut << sv.str() << std::endl;




    //int numPixels = samplePixels.rows; //may differ from sampleSize

    //Get column mean values for the samplePixels
    //cv::Mat columnMeans(1, 3, cv::DataType<double>::type);
    //for (int i = 0; i < 3; i++) {
    //    columnMeans.at<double>(0, i);
    //    
    //sv << cv::mean(samplePixels) << std::endl; //returns a cv::Scalar. now what do I do with that?
    //}


    //Time the operation
    //Start the steady clock
    auto startTime = std::chrono::steady_clock::now();






    //use PCA for now
    cv::PCA pcAnalysis(samplePixels, cv::noArray(), cv::PCA::DATA_AS_ROW, 2); //2 components to project onto a plane
    //Project all of the sampled pixels into the new basis
    cv::Mat projectedPoints;
    pcAnalysis.project(samplePixels, projectedPoints);

    std::stringstream scov;
    scov << "The PCA eigenvalues are: " << pcAnalysis.eigenvalues << std::endl;
    scov << "The PCA basis vectors are: " << pcAnalysis.eigenvectors << std::endl;
    scov << "The first projected point is: " << projectedPoints.row(0) << std::endl;

    //Is projecting into the new basis all that difficult?
    cv::Mat basisVectors;
    cv::transpose(pcAnalysis.eigenvectors, basisVectors);

    scov << "The basisVectors (original directions): " << std::endl;
    scov << basisVectors << std::endl;



    //Here's where the MacenkoHistogram will go. Pass it the points and the unaltered basis vectors
    std::unique_ptr<MacenkoHistogram> histMaker = std::make_unique<MacenkoHistogram>();
    //What should the signs of the basis vectors be? 
    cv::Mat optSignBasisVectors;
    histMaker->OptimizeBasisVectorSigns(samplePixels, basisVectors, optSignBasisVectors, 
        MacenkoHistogram::VectorDirection::COLUMNVECTORS);




    if (basisVectors.at<double>(0, 0) < 0.0) {
        basisVectors.at<double>(0, 0) *= -1.0;
        basisVectors.at<double>(1, 0) *= -1.0;
        basisVectors.at<double>(2, 0) *= -1.0;
    }
    if (basisVectors.at<double>(0, 1) < 0.0) {
        basisVectors.at<double>(0, 1) *= -1.0;
        basisVectors.at<double>(1, 1) *= -1.0;
        basisVectors.at<double>(2, 1) *= -1.0;
    }

    scov << "The basisVectors (new directions): " << std::endl;
    scov << basisVectors << std::endl;

    cv::Mat projectedPointsByMM;
    cv::gemm(samplePixels, basisVectors, 1, cv::Mat(), 0, projectedPointsByMM, 0);

    scov << "The gemm output: " << std::endl;
    scov << projectedPointsByMM << std::endl;

    //Get angular coordinates with respect to the eigenvectors (the new orthogonal basis)
    cv::Mat angularCoords(projectedPointsByMM.rows, 1, cv::DataType<float>::type);

    for (auto p = angularCoords.begin<float>(); p != angularCoords.end<float>(); p++) {
        int row = static_cast<int>(p.lpos());
        bool angleUndef = (projectedPointsByMM.at<float>(row, 0) == 0.0)
            && (projectedPointsByMM.at<float>(row, 1) == 0.0);
        //Use maximum float value as undefined value
        float undefined = std::numeric_limits<float>::max();
        *p = angleUndef ? undefined
            : std::atan2(projectedPointsByMM.at<float>(row, 1), projectedPointsByMM.at<float>(row, 0));
    }

    //Create a histogram of values between -pi and pi
    cv::Mat angleHist;
    int channels[1] = { 0 };
    int histSize[1] = { 1024 };
    float pi = static_cast<float>(CV_PI);
    float range[] = { -pi, pi };
    const float* histRange[] = { range };
    bool uniform = true;
    bool accumulate = false;
    cv::calcHist(&angularCoords, 1, channels, cv::Mat(), angleHist, 1, histSize, histRange, uniform, accumulate);
    //The resulting type of angleHist elements is float

    scov << "This is angleHist size: " << std::endl;
    scov << (angleHist.size)[0] << std::endl;
    scov << angleHist.type() << std::endl;

    //Count how many elements were added to the histogram
    cv::Scalar scalarTotal = cv::sum(angleHist);
    double histoCountTotal = scalarTotal[0];

    //Calculate the slope and intercept to convert from bin to angle
    //angle = intercept + slope*bin
    double intercept = static_cast<double>(range[0]);
    double slope = static_cast<double>(range[1] - range[0]) / static_cast<double>(histSize[0]);

    //Find the bins containing the lower and upper percentiles
    double lowerFraction = (percentileThreshold) / 100.0;
    double upperFraction = (100.0 - percentileThreshold) / 100.0;
    bool lowerPassed(false), upperPassed(false);
    double lowerBin(-1.0), upperBin(-1.0);
    double cumulativeSum(0.0);
    for (int bin = 0; bin < (angleHist.size)[0]; bin++) {
        double prevFraction = cumulativeSum / histoCountTotal;
        cumulativeSum += static_cast<double>(angleHist.at<float>(bin, 0));
        double currentFraction = cumulativeSum / histoCountTotal;

        //scov << bin << " currentFraction: " << currentFraction << std::endl;

        if (!lowerPassed && (currentFraction >= lowerFraction)) {
            lowerPassed = true;
            //linearly interpolate a bin position
            lowerBin = bin + (lowerFraction - prevFraction) / (currentFraction - prevFraction);
        }
        if (!upperPassed && (currentFraction >= upperFraction)) {
            upperPassed = true;
            //linearly interpolate a bin position
            upperBin = bin + (upperFraction - prevFraction) / (currentFraction - prevFraction);
        }
        if (lowerPassed && upperPassed) { break; }
    }
    //Convert from bin values to angles
    double lowerAngle = intercept + slope * lowerBin;
    double upperAngle = intercept + slope * upperBin;


    scov << "Lower bin number: " << lowerBin << ", upper bin number: " << upperBin << std::endl;
    scov << "Lower angle: " << lowerAngle << ", upper angle: " << upperAngle << std::endl;

    //Convert to Cartesian on the PCA projected plane
    cv::Mat stainVectorProjections(2, 2, cv::DataType<double>::type);
    stainVectorProjections.at<double>(0, 0) = std::cos(lowerAngle);
    stainVectorProjections.at<double>(0, 1) = std::sin(lowerAngle);
    stainVectorProjections.at<double>(1, 0) = std::cos(upperAngle);
    stainVectorProjections.at<double>(1, 1) = std::sin(upperAngle);

    //Back-project to get un-normalized stain vectors
    cv::Mat stainVectorOutput;
    pcAnalysis.backProject(stainVectorProjections, stainVectorOutput);


    scov << "Well, here are what it calculates as the stain vectors:" << std::endl;
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

    StainVectorMath::Make3x3MatrixUnitary(nonUnitary, outputVectors);
    scov << "Here they are unitary: " << std::endl;

    scov << outputVectors[0] << ", " << outputVectors[1] << ", " << outputVectors[2] << std::endl;
    scov << outputVectors[3] << ", " << outputVectors[4] << ", " << outputVectors[5] << std::endl;
    scov << outputVectors[6] << ", " << outputVectors[7] << ", " << outputVectors[8] << std::endl;


    tempOut << scov.str() << std::endl;





    //How long did it take to get the basis vectors?
    auto afterBasisVectorsTime = std::chrono::steady_clock::now();

    //How long did getting the basis vectors take?
    long long gettingBasisVectorsDuration = std::chrono::duration_cast<std::chrono::microseconds>(afterBasisVectorsTime - startTime).count();














    //Temp file output
    //std::fstream tempOut;
    //tempOut.open("D:\\mschumaker\\projects\\Sedeen\\testData\\output\\tempout.txt", std::fstream::out);



    //Let's write these to a file to check that something's happening
    //std::stringstream sv;
    //sv << numPixelsAddedToMatrix-1 << ") Tile " << tl << ", Pixel " << px << ": ";
    //sv << sampledPixelsMatrix.at<double>(numPixelsAddedToMatrix-1, 0) << ", ";
    //sv << sampledPixelsMatrix.at<double>(numPixelsAddedToMatrix-1, 1) << ", ";
    //sv << sampledPixelsMatrix.at<double>(numPixelsAddedToMatrix-1, 2);
    //tempOut << sv.str() << std::endl;


    //tempOut << sv.str() << std::endl;
    tempOut.close();




}//end single-parameter ComputeStainVectors






//This overload does not have a default value for sampleSize, so it requires at least two arguments.
//Thus there is a clear difference in arguments between this and the other overload of the method
void StainVectorMacenko::ComputeStainVectors(double outputVectors[9], int sampleSize,
    double ODthreshold /* = 0.15 */, double percentileThreshold /* = 1.0 */) {
    assert(nullptr != this->GetSourceFactory());
    //Set member variables with the argument values
    this->SetSampleSize(sampleSize);
    this->SetODThreshold(ODthreshold);
    this->SetPercentileThreshold(percentileThreshold);
    //Call the single-parameter version of this method, which uses the member variables
    this->ComputeStainVectors(outputVectors);
}//end multi-parameter ComputeStainVectors

} // namespace image
} // namespace sedeen
