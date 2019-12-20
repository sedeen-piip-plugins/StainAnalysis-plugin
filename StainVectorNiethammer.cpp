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

#include "StainVectorNiethammer.h"

//#include <chrono>
//#include <random>
#include <fstream>
#include <sstream>

//For now, include ODConversion here, but try to do the OD conversion and thresholding
//in a kernel, and use a factory to apply it before passing the factory to this class
#include "ODConversion.h"
#include "StainVectorMath.h"
#include "NiethammerHistogram.h"
#include "BasisTransform.h"

namespace sedeen {
namespace image {

StainVectorNiethammer::StainVectorNiethammer(std::shared_ptr<tile::Factory> source,
    double ODthreshold /*= 0.15 */, double percentileThreshold /*= 1.0 */, double qAdjustmentFactor /*= 0.15 */)
    : StainVectorOpenCV(source),
    m_sampleSize(0), //Must set to greater than 0 to ComputeStainVectors
    m_avgODThreshold(ODthreshold),     //assign default value
    m_percentileThreshold(percentileThreshold), //assign default value
    m_qVectorMixRatio(qAdjustmentFactor),    //assign default value
    m_priors({ 0.0 })           //prior values 0 unless redefined after construction
{}//end constructor

StainVectorNiethammer::~StainVectorNiethammer(void) {
}//end destructor

void StainVectorNiethammer::ComputeStainVectors(double (&outputVectors)[9]) {
    if (this->GetSourceFactory() == nullptr) { return; }
    //Using this overload of the method requires the sample size to be already set
    int sampleSize = this->GetSampleSize();
    if (sampleSize <= 0) { return; }
    //This overload requires the ODthreshold to be already set
    double ODthreshold = this->GetODThreshold();
    //This overload requires the percentileThreshold to be already set
    double percentileThreshold = this->GetPercentileThreshold();
    if (percentileThreshold <= 0.0) { return; }
    //This overload requires that the q adjustment factor to be already set
    double qVectorMixRatio = this->GetQVectorMixRatio();
    //Priors may be zero, or must be set already for this method overload
    double thePriors[9];
    this->GetPriors(thePriors);
    //Convert priors to CV Mat
    cv::Mat cvPriors;
    this->StainCArrayToCVMat(thePriors, cvPriors, true, 2); //number of rows to output


    //Sample a set of pixel values from the source
    cv::Mat samplePixels;
    auto theSampler = this->GetRandomWSISampler();
    if (theSampler == nullptr) { return; }
    bool samplingSuccess = theSampler->ChooseRandomPixels(samplePixels, sampleSize, ODthreshold);
    if (!samplingSuccess) { return; }

    //Temp file output
    std::fstream tempOut;
    tempOut.open("D:\\mschumaker\\projects\\Sedeen\\testData\\output\\tempout-ComputeNiethammer.txt", std::fstream::out);
    std::stringstream scov;

    //Create a mask to identify cluster assignments, and a holder for the values from the previous iteration
    //0 and 1 are cluster assignments. Initialize to cluster 0.
    cv::Size sizeOfMask = cv::Size(samplePixels.rows, 1);
    int ctype = cv::DataType<int>::type;
    cv::Mat clusterAssignments = cv::Mat::zeros(sizeOfMask, ctype);
    cv::Mat prevClusterAssignments; //initially empty

    //I guess do stuff here and then refactor when I know what I'm doing?

    //Steps:
    //project onto basis defined by the priors, go ahead and find the mean, though.
    std::unique_ptr<BasisTransform> testingBasisTransform = std::make_unique<BasisTransform>();
    cv::Mat projectedPoints;
    testingBasisTransform->NiethammerProjection(samplePixels, projectedPoints, cvPriors);
    
    //Compute q1 and q2 by mixing the stain priors
    cv::Mat qVectors, projQPriors;
    ComputeQVectorsFromPriors(cvPriors, qVectors, this->GetQVectorMixRatio());

    //TEMP: create a separate basis transform object for projection of the q priors
    std::unique_ptr<BasisTransform> qBasisTransform = std::make_unique<BasisTransform>();
    qBasisTransform->NiethammerProjection(qVectors, projQPriors, cvPriors);

    scov << "The original cvPriors: " << cvPriors << std::endl;
    scov << "The adjusted qPriors:  " << qVectors << std::endl;
    scov << "The projected qPriors: " << projQPriors << std::endl;


    //Create a histogramming class that can identify clusters
    std::unique_ptr<NiethammerHistogram> testingHistogram = std::make_unique<NiethammerHistogram>();

    //Assign clusterAssignments reference to prevClusterAssignments, get new clusterAssignments
    prevClusterAssignments = clusterAssignments;
    testingHistogram->AssignClusters(projectedPoints, clusterAssignments, qVectors);
    

    //Test equality of the previous and the new cluster assignment matrices
    bool assignmentsEqual = AreEqual(prevClusterAssignments, clusterAssignments);

    scov << "Are the old and new cluster assignments equivalent? " << assignmentsEqual << std::endl;


    ////Create a class to perform the basis transformation of the sample pixels
    //std::unique_ptr<BasisTransform> theBasisTransform = std::make_unique<BasisTransform>();
    ////Both the input and output data points should be the matrix rows (columns are pixel elements)
    //cv::Mat projectedPoints;
    //theBasisTransform->PCAPointTransform(samplePixels, projectedPoints);

    //Now, the only reason to see the basis vectors in this class is out of curiosity
    //cv::Mat basisVectors;
    //bool getBasisSuccess = theBasisTransform->GetBasisVectors(basisVectors);

    //scov << "The basis vectors are: " << basisVectors << std::endl;
    //scov << "the projectedPoints are: " << std::endl;
    //scov << projectedPoints << std::endl;




    ////Create a class to histogram the results and find 2D vectors corresponding to percentile thresholds
    //std::unique_ptr<MacenkoHistogram> theHistogram = std::make_unique<MacenkoHistogram>();
    //cv::Mat percentileThreshVectors;
    //theHistogram->PercentileThresholdVectors(projectedPoints, percentileThreshVectors, this->GetPercentileThreshold());

    //scov << "The percentile threshold vectors: " << std::endl;
    //scov << percentileThreshVectors << std::endl;


    ////Back-project to get un-normalized stain vectors. DO NOT translate to the mean after backprojection.
    //cv::Mat backProjectedVectors;
    //theBasisTransform->backProjectPoints(percentileThreshVectors, backProjectedVectors, false); //useMean=false

    //scov << "The back projected percentileThreshVectors: " << std::endl;
    //scov << backProjectedVectors << std::endl;


    ////Convert to C array and normalize rows
    //double tempStainVecOutput[9] = {0.0};
    //StainCVMatToCArray(backProjectedVectors, tempStainVecOutput, true);
    //std::copy(std::begin(tempStainVecOutput), std::end(tempStainVecOutput), std::begin(outputVectors));




    //scov << "Testing CVMat to C array conversion (normalize=true): " << std::endl;
    //for (int i = 0; i < 9; i++) {
    //    scov << tempStainVecOutput[i] << ", ";
    //}
    //scov << std::endl;
    



    //Test converting it the other way
    //cv::Mat convertBack;
    //StainCArrayToCVMat(tempStainVecOutput, convertBack, true);

    //scov << "Testing C array to CVMat conversion (normalize=true): " << std::endl;
    //scov << convertBack << std::endl;


    tempOut << scov.str() << std::endl;
    tempOut.close();


}//end single-parameter ComputeStainVectors





//This overload does not have a default value for sampleSize, so it requires two arguments
void StainVectorNiethammer::ComputeStainVectors(double(&outputVectors)[9], const int sampleSize) {
    //This method creates a priors array initialized to zeros, and calls the method with priors
    double zeroPriors[9] = { 0.0 };
    this->ComputeStainVectors(outputVectors, zeroPriors, sampleSize);
}//end no-priors multi-parameter ComputeStainVectors

//This overload requires three arguments
void StainVectorNiethammer::ComputeStainVectors(double (&outputVectors)[9], 
    const double (&inputPriors)[9], const int sampleSize) {
    if (this->GetSourceFactory() == nullptr) { return; }
    //Set member variables with the argument values
    this->SetSampleSize(sampleSize);
    //Set the prior stain vector values
    this->SetPriors(inputPriors);
    //Call the single-parameter version of this method, which uses the member variables
    this->ComputeStainVectors(outputVectors);
}//end multi-parameter ComputeStainVectors

void StainVectorNiethammer::ComputeQVectorsFromPriors(cv::InputArray stainPriors, cv::OutputArray qVectors, double qAdjustmentFactor) {
    //The q adjustment factor should be between 0 and 0.5, but that will not be enforced
    if (stainPriors.empty()) { return; }
    int inputRows = stainPriors.rows();
    int inputCols = stainPriors.cols();
    int ctype = stainPriors.depth();
    cv::Mat stainPriorMat = stainPriors.getMat();
    cv::Mat qVectorMat(inputRows, inputCols, ctype);

    //This method assumes row vectors
    //If there is only one row, assign input to output, return
    if (inputRows == 1) {
        qVectorMat = stainPriorMat.clone();
        qVectors.assign(qVectorMat);
    }
    else {
        //Get pointers to the first two rows (though there should only be two)
        cv::Mat s1 = stainPriorMat.row(0);
        cv::Mat s2 = stainPriorMat.row(1);
        //Get pointers to the destination rows
        cv::Mat q1 = qVectorMat.row(0);
        cv::Mat q2 = qVectorMat.row(1);
        //alpha = qAdjustmentFactor
        //q1 = (1-alpha)s1 + (alpha)s2
        //q2 = (alpha)s1 + (1-alpha)s2
        q1 = s1.mul(1.0 - qAdjustmentFactor) + s2.mul(qAdjustmentFactor);
        q2 = s1.mul(qAdjustmentFactor) + s2.mul(1.0 - qAdjustmentFactor);
        qVectors.assign(qVectorMat);
    }
}//end ComputeQVectorsFromPriors

} // namespace image
} // namespace sedeen
