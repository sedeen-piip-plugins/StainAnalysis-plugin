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

StainVectorMacenko::StainVectorMacenko(std::shared_ptr<tile::Factory> source,
    double ODthreshold /* = 0.15 */, double percentileThreshold /* = 1.0 */,
    int numHistoBins /* = 1024 */)
    : StainVectorOpenCV(source),
    m_sampleSize(0), //Must set to greater than 0 to ComputeStainVectors
    m_avgODThreshold(ODthreshold), //assign default value
    m_percentileThreshold(percentileThreshold), //assign default value
    m_numHistogramBins(numHistoBins) //assign default value
{}//end constructor

StainVectorMacenko::~StainVectorMacenko(void) {
}//end destructor

void StainVectorMacenko::ComputeStainVectors(double (&outputVectors)[9]) {
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
    std::unique_ptr<MacenkoHistogram> theHistogram
        = std::make_unique<MacenkoHistogram>(this->GetPercentileThreshold(), this->GetNumHistogramBins());
    cv::Mat percentileThreshVectors;
    theHistogram->PercentileThresholdVectors(projectedPoints, percentileThreshVectors);

    scov << "The percentile threshold vectors: " << std::endl;
    scov << percentileThreshVectors << std::endl;


    //Back-project to get un-normalized stain vectors. DO NOT translate to the mean after backprojection.
    cv::Mat backProjectedVectors;
    theBasisTransform->backProjectPoints(percentileThreshVectors, backProjectedVectors, false); //useMean=false

    scov << "The back projected percentileThreshVectors: " << std::endl;
    scov << backProjectedVectors << std::endl;


    //Convert to C array and normalize rows
    double tempStainVecOutput[9] = {0.0};
    StainCVMatToCArray(backProjectedVectors, tempStainVecOutput, true);
    std::copy(std::begin(tempStainVecOutput), std::end(tempStainVecOutput), std::begin(outputVectors));




    scov << "Testing CVMat to C array conversion (normalize=true): " << std::endl;
    for (int i = 0; i < 9; i++) {
        scov << tempStainVecOutput[i] << ", ";
    }
    scov << std::endl;
    



    //Test converting it the other way
    //cv::Mat convertBack;
    //StainCArrayToCVMat(tempStainVecOutput, convertBack, true);

    //scov << "Testing C array to CVMat conversion (normalize=true): " << std::endl;
    //scov << convertBack << std::endl;


    tempOut << scov.str() << std::endl;
    tempOut.close();





}//end single-parameter ComputeStainVectors



//This overload does not have a default value for sampleSize, so it requires at two arguments
void StainVectorMacenko::ComputeStainVectors(double (&outputVectors)[9], const int sampleSize) {
    if (this->GetSourceFactory() == nullptr) { return; }
    //Set member variables with the argument values
    this->SetSampleSize(sampleSize);
    //Call the single-parameter version of this method, which uses member variables
    this->ComputeStainVectors(outputVectors);
}//end multi-parameter ComputeStainVectors

} // namespace image
} // namespace sedeen
