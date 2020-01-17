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

#include "StainVectorICA.h"

//The RADICAL algorithm for independent component analysis (ICA)
#include <mlpack/methods/radical/radical.hpp>

//#include <chrono>
//#include <random>
#include <fstream>
#include <sstream>

//For now, include ODConversion here, but try to do the OD conversion and thresholding
//in a kernel, and use a factory to apply it before passing the factory to this class
#include "ODConversion.h"
#include "StainVectorMath.h"

namespace sedeen {
namespace image {

StainVectorICA::StainVectorICA(std::shared_ptr<tile::Factory> source,
    double ODthreshold /*= 0.15 */)
    : StainVectorMLPACK(source),
    m_sampleSize(0), //Must set to greater than 0 to ComputeStainVectors
    m_numStains(2),  //Can be 2 or 3
    m_avgODThreshold(ODthreshold) //assign default value
{}//end constructor

StainVectorICA::~StainVectorICA(void) {
}//end destructor

void StainVectorICA::ComputeStainVectors(double (&outputVectors)[9]) {
    if (this->GetSourceFactory() == nullptr) { return; }
    //Using this overload of the method requires setting sample size in advance
    int sampleSize = this->GetSampleSize();
    if (sampleSize <= 0) { return; }
    double ODthreshold = this->GetODThreshold();

    //Sample a set of pixel values from the source
    cv::Mat samplePixels;
    auto theSampler = this->GetRandomWSISampler();
    if (theSampler == nullptr) { return; }
    bool samplingSuccess = theSampler->ChooseRandomPixels(samplePixels, sampleSize, ODthreshold);
    if (!samplingSuccess) { return; }

    //Convert samplePixels from CV to Armadillo
    arma::Mat<double> armaSamplePixels = CVMatToArmaMat<double>(samplePixels);
    arma::Mat<double> transposedArmaPixels = arma::trans(armaSamplePixels);


    //Implement the RADICAL algorithm
    //RADICAL has a number of parameters with default values
    //Define them here for clarity and potential tuning
    //Descriptions from https://rcppmlpack.github.io/mlpack-doxygen/classmlpack_1_1radical_1_1Radical.html
    const double noiseStdDev = 0.175; //Std deviation of the Gaussian noise added to the replicates of the data points during Radical2D
    const size_t replicates = 30;     //Number of Gaussian-perturbed replicates to use (per point) in Radical2D
    const size_t angles = 150;        //Number of angles to consider in brute-force search during Radical2D
    const size_t sweeps = 0;          //Number of sweeps. Each sweep calls Radical2D once for each pair of dimensions
    const size_t m = 0;               //The variable m from Vasicek's m-spacing estimator of entropy.
    mlpack::radical::Radical radicalICA(noiseStdDev, replicates, angles, sweeps, m);
    arma::Mat<double> independentComponents, unmixingMatrix;
    radicalICA.DoRadical(transposedArmaPixels, independentComponents, unmixingMatrix);


    //Ok, I can run ICA, but the results aren't useful to me. Can I improve on this?



    //The stain values are in the encoding matrix. Convert to output array
    //cv::Mat encodingAsCV = ArmaMatToCVMat<double>(encodingMat);


    //Temp file output
    std::fstream tempOut;
    tempOut.open("D:\\mschumaker\\projects\\Sedeen\\testData\\output\\tempout-ComputeStainVectorsICA.txt", std::fstream::out);

    std::stringstream scov;

    scov << "These are the inputs in CV format: " << std::endl;
    scov << samplePixels << std::endl;
    scov << "\n\nThese are the transposed inputs in Armadillo format: " << std::endl;
    scov << transposedArmaPixels << std::endl;


    scov << "\n\nSo what is the ICA output? " << std::endl;
    scov << "independentComponents: " << independentComponents << std::endl;
    scov << "unmixingMatrix: " << unmixingMatrix << std::endl;

    //Convert to C array and normalize rows
    double tempStainVecOutput[9] = {0.0};
    //StainCVMatToCArray(unmixingAsCV, tempStainVecOutput, true);
    for (int i = 0; i < 9; i++) {
        outputVectors[i] = tempStainVecOutput[i];
    }


    tempOut << scov.str() << std::endl;
    tempOut.close();

}//end single-parameter ComputeStainVectors



//This overload does not have a default value for sampleSize, so it requires at two arguments
void StainVectorICA::ComputeStainVectors(double (&outputVectors)[9], int sampleSize) {
    if (this->GetSourceFactory() == nullptr) { return; }
    //Set member variables with the argument values
    this->SetSampleSize(sampleSize);
    //Call the single-parameter version of this method, which uses the member variables
    this->ComputeStainVectors(outputVectors);
}//end multi-parameter ComputeStainVectors

} // namespace image
} // namespace sedeen
