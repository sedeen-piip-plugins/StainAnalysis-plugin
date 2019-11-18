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

#include "BasisTransform.h"

#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>

namespace sedeen {
namespace image {

BasisTransform::BasisTransform() : m_numTestingPixels(10),
    m_rgen((std::random_device())()) //Initialize random number generation
{
}//end constructor

BasisTransform::~BasisTransform(void) {
}//end destructor

void BasisTransform::PCAPointTransform(cv::InputArray sourcePoints, cv::OutputArray outputPoints, 
    cv::InputArray sourceMask /*=cv::noArray()*/, cv::InputArray inputMean /*=cv::noArray()*/, 
    const VectorDirection &sourcePointDir /*= VectorDirection::ROWVECTORS*/) {
    if (sourcePoints.empty()) { return; }

    //Get the source points as a matrix
    cv::Mat sourceMat(sourcePoints.getMat());
    //Get the input mean as a matrix
    cv::Mat inputMeanMat(inputMean.getMat());

    //Determine the data type of the source data
    int ctype = sourceMat.depth();

    int covar_flags = cv::CovarFlags::COVAR_SCALE | cv::CovarFlags::COVAR_NORMAL;
    int numElements, numPoints;
    cv::Size sizeOfMean;

    //If sourcePointDir == VectorDirection::ROWVECTORS, data points are entered as rows, RGB elements columns
    if (sourcePointDir == VectorDirection::ROWVECTORS) {
        numElements = sourceMat.cols;
        numPoints = sourceMat.rows;
        sizeOfMean = cv::Size(numElements, 1);
        covar_flags |= cv::CovarFlags::COVAR_ROWS;
    }
    else { //VectorDirection::COLUMNVECTORS
        numElements = sourceMat.rows;
        numPoints = sourceMat.cols;
        sizeOfMean = cv::Size(1, numElements);
        covar_flags |= cv::CovarFlags::COVAR_COLS;
    }

    //We will only consider over-determined cases in this class: numPoints > numElements
    if (numPoints <= numElements) { return; }

    //Define matrix for element means
    cv::Mat elementMeans(sizeOfMean, ctype);
    //Define and allocate space for covariance matrix
    cv::Mat covar(numElements, numElements, ctype);

    //If inputMeanMat is not empty, copy to elementMeans (if same size)
    if (!inputMeanMat.empty()) {
        //If the size does not match elementMeans, stop and return
        if (inputMeanMat.size() != sizeOfMean) { return; }
        inputMeanMat.convertTo(elementMeans, ctype); //converts and copies to elementMeans
        //Set the bit to use pre-calculated means in the covariance matrix
        covar_flags |= cv::CovarFlags::COVAR_USE_AVG;
    }

    //Calculate the covariance matrix
    cv::calcCovarMatrix(sourceMat, covar, elementMeans, covar_flags, ctype);

    //Calculate the eigenvalues and eigenvectors
    cv::Mat eigenvalues, eigenvectors;
    cv::eigen(covar, eigenvalues, eigenvectors);

    //Set the mean, eigenvalues, eigenvectors, and basis vectors using covar and eigen outputs
    SetPointMean(elementMeans);
    SetEigenvalues(eigenvalues);
    SetEigenvectors(eigenvectors);


    //Here is where to optimize the basis vector directions, if we still want to do that
    //cv::Mat optSignBasisVectors;
    //OptimizeBasisVectorSigns(sourcePoints, basisVectors, optSignBasisVectors,
    //    BasisTransform::VectorDirection::COLUMNVECTORS);
    //???
    SetBasisVectors(GetEigenvectors(m_reqdBasisVectors, GetEigenvectorElementsDirection()), GetEigenvectorElementsDirection());


    std::fstream tempOut;
    tempOut.open("D:\\mschumaker\\projects\\Sedeen\\testData\\output\\tempout-PCAPointTransform.txt", std::fstream::out);
    std::stringstream sv;
    sv << "Am I getting the projection right?" << std::endl;

    sv << "These are the raw source points: " << std::endl;
    sv << sourceMat << std::endl;

    sv << "This is the covariance matrix from calcCovarMatrix: " << std::endl;
    sv << covar << std::endl;

    sv << "The eigen vectors from covar and eigen are: " << std::endl;
    sv << this->GetEigenvectors() << std::endl;

    sv << "The basis vectors from covar and eigen are: " << std::endl;
    sv << this->GetBasisVectors() << std::endl;

    sv << "This is the mean values covar is coming up with: " << std::endl;
    sv << this->GetPointMean() << std::endl;

    //Project to put vectors on a plane. DO NOT translate to the mean before projection.
    cv::Mat projectedPointsOut;
    bool useMean = false;
    this->projectPoints(sourcePoints, projectedPointsOut, useMean);
    sv << "The projected points from my projectPoints method are: " << std::endl;
    sv << projectedPointsOut << std::endl;

    sv << "************************************************" << std::endl;


    //cv::Mat backProjTest;
    //bool backProjSuccess = this->backProjectPoints(projectedPointsOut, backProjTest, useMean);
    //sv << "These are back projected points from my method: " << std::endl;
    //sv << backProjTest << std::endl;
    //sv << "This is the value of backProjSuccess: " << backProjSuccess << std::endl;






    tempOut << sv.str() << std::endl;
    tempOut.close();

    //Assign projectedPoints to the output
    outputPoints.assign(projectedPointsOut);
}//end PCAPointTransform




bool BasisTransform::projectPoints(cv::InputArray sourcePoints, cv::OutputArray projectedPoints, bool subtractMean /*= true*/) const {
    cv::Mat basisVecs, means, tempProjPoints;
    //Need the basis vectors and point element means
    this->GetPointMean(means);
    bool getVecsSuccess = this->GetBasisVectors(basisVecs);
    if (!getVecsSuccess || means.empty()) { return false; }

    this->projectPoints(sourcePoints, tempProjPoints, basisVecs, means, subtractMean);
    if (tempProjPoints.empty()) {
        return false;
    }
    else {
        projectedPoints.assign(tempProjPoints);
        return true;
    }
    return false;
}//end projectPoints, public version

void BasisTransform::projectPoints(cv::InputArray sourcePoints, cv::OutputArray projectedPoints, 
    cv::InputArray basisVectors, cv::InputArray means, bool subtractMean /*= true*/) const {
    //Check that the sourcePoints, basisVectors matrices aren't empty
    if (sourcePoints.empty() || basisVectors.empty()) { return; }
    //Even if subtractMean is false, a means matrix is still required to deduce matrix orientation
    if (means.empty()) { return; }
    //Check that basisVectors isn't all zeros
    if (cv::countNonZero(basisVectors) == 0) { return; }

    //Check that matrix dimensions are compatible
    bool validSizes = (means.rows() == 1 && means.cols() == sourcePoints.cols())
        || (means.cols() == 1 && means.rows() == sourcePoints.rows());
    if (!validSizes) { return; }

    cv::Mat sourceMat(sourcePoints.getMat());
    int sourceType = sourceMat.type();
    //Match the types of the basis vectors and means to the source point type
    cv::Mat meansMat, _meansMat(means.getMat());
    _meansMat.convertTo(meansMat, sourceType);
    cv::Mat basisMat, _basisMat(basisVectors.getMat());
    _basisMat.convertTo(basisMat, sourceType);

    //If subtractMean is true, fill a matrix with the mean values to subtract from data elements
    //If subtractMean is false, fill the matrix with zeros
    int yReps = sourceMat.rows / meansMat.rows;
    int xReps = sourceMat.cols / meansMat.cols;
    cv::Mat repeatedMeans;
    if (subtractMean) {
        repeatedMeans = cv::repeat(meansMat, yReps, xReps);
    }
    else {
        cv::Mat zeroVals = cv::Mat::zeros(meansMat.rows, meansMat.cols, meansMat.type());
        repeatedMeans = cv::repeat(zeroVals, yReps, xReps);
    }

    //Subtract the means from the source data matrix
    cv::Mat sourceMinusMeans;
    cv::subtract(sourceMat, repeatedMeans, sourceMinusMeans);

    //Determine the multiplication order by the orientation of the mean matrix
    if (meansMat.rows == 1) {
        cv::gemm(sourceMinusMeans, basisMat, 1, cv::Mat(), 0, projectedPoints, cv::GemmFlags::GEMM_2_T);
    }
    else if (meansMat.cols == 1) {
        cv::gemm(basisMat, sourceMinusMeans, 1, cv::Mat(), 0, projectedPoints, 0);
    }
}//end projectPoints, protected 4-argument version

bool BasisTransform::backProjectPoints(cv::InputArray projectedPoints, cv::OutputArray backProjPoints, bool addMean /*= true*/) const {
    cv::Mat basisVecs, means, tempBackProjPoints;
    //Need the basis vectors and point element means
    this->GetPointMean(means);
    //If addMean is true, check that the means matrix is not empty
    if (addMean && means.empty()) { return false; }
    //Check that basis vectors are not empty
    bool getVecsSuccess = this->GetBasisVectors(basisVecs);
    if (!getVecsSuccess || means.empty()) { return false; }

    this->backProjectPoints(projectedPoints, tempBackProjPoints, basisVecs, means);
    if (tempBackProjPoints.empty()) {
        return false;
    }
    else {
        backProjPoints.assign(tempBackProjPoints);
        return true;
    }
    return false;
}//end backProjectPoints, public version

void BasisTransform::backProjectPoints(cv::InputArray projectedPoints, cv::OutputArray backProjPoints, 
    cv::InputArray basisVectors, cv::InputArray means, bool addMean /*= true*/) const {
    //Check that the projectedPoints, basisVectors matrices aren't empty
    if (projectedPoints.empty() || basisVectors.empty()) { return; }
    //Even if addMean is false, a means matrix is still required to deduce matrix orientation
    if (means.empty()) { return; }
    //Check that basisVectors isn't all zeros
    if (cv::countNonZero(basisVectors) == 0) { return; }

    //Check that matrix dimensions are compatible
    bool validSizes = (means.rows() == 1 && basisVectors.rows() == projectedPoints.cols())
        || (means.cols() == 1 && basisVectors.rows() == projectedPoints.rows());
    if (!validSizes) { return; }

    cv::Mat projMat(projectedPoints.getMat());
    int projType = projMat.type();
    //Match the types of the basis vectors and means to the data point type
    cv::Mat meansMat, _meansMat(means.getMat());
    _meansMat.convertTo(meansMat, projType);
    cv::Mat basisMat, _basisMat(basisVectors.getMat());
    _basisMat.convertTo(basisMat, projType);

    //If addMean is true, fill a matrix with mean values to add after back projection
    //If addMean is false, fill the matrix with zeros
    if (addMean) {
        //Determine the multiplication order by the orientation of the mean matrix
        if (meansMat.rows == 1) {
            cv::Mat repeatedMeans = cv::repeat(meansMat, projMat.rows, 1);
            cv::gemm(projMat, basisMat, 1, repeatedMeans, 1, backProjPoints, 0);
        }
        else if (meansMat.cols == 1) {
            cv::Mat repeatedMeans = cv::repeat(meansMat, 1, projMat.cols);
            cv::gemm(basisMat, projMat, 1, repeatedMeans, 1, backProjPoints, cv::GemmFlags::GEMM_1_T);
        }
    }
    else {
        cv::Mat zeroVals = cv::Mat::zeros(meansMat.rows, meansMat.cols, meansMat.type());
        //Determine the multiplication order by the orientation of the mean matrix
        if (meansMat.rows == 1) {
            cv::Mat repeatedMeans = cv::repeat(zeroVals, projMat.rows, 1);
            cv::gemm(projMat, basisMat, 1, repeatedMeans, 1, backProjPoints, 0);
        }
        else if (meansMat.cols == 1) {
            cv::Mat repeatedMeans = cv::repeat(zeroVals, 1, projMat.cols);
            cv::gemm(basisMat, projMat, 1, repeatedMeans, 1, backProjPoints, cv::GemmFlags::GEMM_1_T);
        }
    }
}//end backProjectPoints, protected 4-argument version

void BasisTransform::SetBasisVectors(cv::InputArray basisVectors, 
    const VectorDirection &evecDir /*= VectorDirection::ROWVECTORS*/) {
    cv::Mat _bvecs = basisVectors.getMat();
    cv::Mat bVecs;
    if (!_bvecs.empty()) {
        if (evecDir == VectorDirection::COLUMNVECTORS) {
            cv::transpose(_bvecs, bVecs);
        }
        else if (evecDir == VectorDirection::ROWVECTORS) {
            bVecs = _bvecs.clone();
        }
        else {
            //Set basis vectors as given
            bVecs = _bvecs.clone();
        }
    }
    this->m_basisVectors = bVecs;
}//end SetBasisVectors

bool BasisTransform::GetBasisVectors(cv::OutputArray basisVectors) const {
    if (this->m_basisVectors.empty()) {
        return false;
    }
    //else
    basisVectors.assign(this->m_basisVectors);
    return true;
}//end GetBasisVectors

cv::Mat BasisTransform::GetBasisVectors() const {
    cv::Mat bVecs;
    bool success = GetBasisVectors(bVecs);
    return bVecs;
}//end GetBasisVectors

void BasisTransform::SetPointMean(cv::InputArray mean) {
    cv::Mat meanMat = mean.getMat();
    if (!meanMat.empty()) {
        this->m_pointMean = meanMat.clone();
    }
}//end SetPointMean

void BasisTransform::GetPointMean(cv::OutputArray mean) const {
    cv::Mat tempMean = m_pointMean.clone();
    mean.assign(tempMean);
}//end GetPointMean

cv::Mat BasisTransform::GetPointMean() const {
    cv::Mat mean;
    GetPointMean(mean);
    return mean;
}//end GetPointMean

void BasisTransform::SetEigenvalues(cv::InputArray evals) {
    cv::Mat evalsMat = evals.getMat();
    if (!evalsMat.empty()) {
        this->m_eigenvalues = evalsMat.clone();
    }
}//end SetEigenvalues

void BasisTransform::GetEigenvalues(cv::OutputArray evals, const int nVals /*= -1*/) const {
    cv::Mat tempVals = m_eigenvalues.clone();
    if (tempVals.empty() || nVals < 0) {
        evals.assign(tempVals);
        return; 
    }
    else {
        //Allow eigenvalues to be either a row or column vector
        int nRows = tempVals.rows;
        int nCols = tempVals.cols;
        cv::Mat subMatrix;
        if (nRows == 1 && nCols > 0) {
            int endCol = (nVals > nCols) ? nCols : nVals;
            subMatrix = tempVals.colRange(0, endCol);
        }
        else if (nRows > 0 && nCols == 1) {
            int endRow = (nVals > nRows) ? nRows : nVals;
            subMatrix = tempVals.rowRange(0, endRow);
        }
        evals.assign(subMatrix);
        return;
    }
}//end GetEigenvalues

cv::Mat BasisTransform::GetEigenvalues(const int nVals /*= -1*/) const {
    cv::Mat evals;
    GetEigenvalues(evals, nVals);
    return evals;
}//end GetEigenvalues

const BasisTransform::VectorDirection BasisTransform::GetEigenvectorElementsDirection() const {
    //Use the eigenvalue axial orientation to determine the orientation of the eigenvectors
    cv::Mat tempVals = m_eigenvalues.clone();
    VectorDirection direction;
    if (tempVals.empty()) {
        direction = VectorDirection::UNDETERMINED;
    }
    else {
        int nRows = tempVals.rows;
        int nCols = tempVals.cols;
        if (nRows == 1 && nCols > 0) {
            direction = VectorDirection::COLUMNVECTORS;
        }
        else if (nRows > 0 && nCols == 1) {
            direction = VectorDirection::ROWVECTORS;
        }
        else {
            direction = VectorDirection::UNDETERMINED;
        }
    }
    return direction;
}//end GetEigenvectorElementsDirection

void BasisTransform::SetEigenvectors(cv::InputArray evecs) {
    cv::Mat evecsMat = evecs.getMat();
    if (!evecsMat.empty()) {
        this->m_eigenvectors = evecsMat.clone();
    }
}//end SetEigenvectors

void BasisTransform::GetEigenvectors(cv::OutputArray evecs, const int nVecs /*= -1*/,
    const VectorDirection &evecDir /*= VectorDirection::ROWVECTORS*/) const {
    //evecDir allows user to specify the vector direction in the Mat
    //The default direction is vectors as rows
    cv::Mat tempVals = m_eigenvectors.clone();
    cv::Mat subMatrix;
    if (tempVals.empty() || nVecs < 0) {
        evecs.assign(tempVals);
        return;
    }
    else {
        int nRows = tempVals.rows;
        int nCols = tempVals.cols;
        if (evecDir == VectorDirection::ROWVECTORS) {
            int endRow = (nVecs > nRows) ? nRows : nVecs;
            subMatrix = tempVals.rowRange(0, endRow);
            evecs.assign(subMatrix);
            return;
        }
        else if (evecDir == VectorDirection::COLUMNVECTORS) {
            int endCol = (nVecs > nCols) ? nCols : nVecs;
            subMatrix = tempVals.colRange(0, endCol);
            evecs.assign(subMatrix);
            return;
        }
        else {
            //Error case: return full matrix
            evecs.assign(tempVals);
            return;
        }
    }
    evecs.assign(tempVals);
}//end GetEigenvectors

cv::Mat BasisTransform::GetEigenvectors(const int nVecs /*= -1*/,
    const VectorDirection &evecDir /*= VectorDirection::ROWVECTORS*/) const {
    cv::Mat evecs;
    GetEigenvectors(evecs, nVecs, evecDir);
    return evecs;
}//end GetEigenvectors





void BasisTransform::OptimizeBasisVectorSigns(cv::InputArray sourcePoints, /*assume sourcePoints to be row vectors */
    cv::InputArray inputVectors, cv::OutputArray outputVectors, bool useMean /*= true*/,
    const VectorDirection &basisVecDir /*= VectorDirection::COLUMNVECTORS*/) {
    //Check a small number of source pixels to try to get projected points with all positive elements
    //The number of pixels to check is arbitrary, default value 10
    cv::Mat subsampleofPixels;


    //pick up here!!!


    //HACK!
    //int numTestPixels = this->GetNumTestingPixels();
    int numTestPixels = 0; 
    
    
    
    if (numTestPixels < 1) {
        //Deep copy the inputVectors to the outputVectors
        outputVectors.assign(inputVectors.getMat().clone());
        return;
    }
    else {
        CreatePixelSubsample(sourcePoints, subsampleofPixels, numTestPixels);
    }

    //Arrange the basis vectors as the columns, if they're not already
    cv::Mat columnBasisVectors;
    if (basisVecDir == VectorDirection::COLUMNVECTORS) {
        //no change
        columnBasisVectors = inputVectors.getMat().clone();
    }
    else if (basisVecDir == VectorDirection::ROWVECTORS) {
        //transpose the inputVectors matrix
        cv::transpose(inputVectors, columnBasisVectors);
    }
    else {
        //invalid, but use the inputVectors as given and proceed anyway
        columnBasisVectors = inputVectors.getMat().clone();
    }

    //Create the list of +/- options for the basis vectors (0 is +, 1 is -)
    int numCombinations = static_cast<int>(std::pow(2, columnBasisVectors.cols));
    cv::Mat testMultCombinations(numCombinations, columnBasisVectors.cols, cv::DataType<int>::type);
    for (int row = 0; row < numCombinations; row++) {
        for (int col = 0; col < columnBasisVectors.cols; col++) {
            int newVal = (row >> col) & 1; //bit shift and mask
            testMultCombinations.at<int>(row, col) = newVal;
        }
    }



    //Temp file output
    std::fstream tempOut;
    tempOut.open("D:\\mschumaker\\projects\\Sedeen\\testData\\output\\tempout-samplepix.txt", std::fstream::out);
    std::stringstream ss;


    //Loop through each of the combinations
    //Get average of subsampled points projected into each variation of the basis vectors
    cv::Mat projAvgPointsByCombo;
    for (int combo = 0; combo < testMultCombinations.rows; combo++) {
        //columnBasisVectors has the basis vectors as columns
        cv::Mat signedBasisVectors = columnBasisVectors.clone();
        for (int vec = 0; vec < testMultCombinations.cols; vec++) {
            double multFactor = (testMultCombinations.at<int>(combo, vec) == 0) ? 1.0 : -1.0;
            signedBasisVectors.col(vec) *=  multFactor;
        }
        ss << "New combo: " << std::endl;
        ss << signedBasisVectors << std::endl;

        //Project the subsample of pixels into this basis

        cv::Mat projectedSamplePoints;
        cv::gemm(subsampleofPixels, signedBasisVectors, 1, cv::Mat(), 0, projectedSamplePoints, 0);

        ss << "The projected subsample of pixels: " << std::endl;
        ss << projectedSamplePoints << std::endl;

        //Use reduce to get column averages
        cv::Mat columnAvg;
        cv::reduce(projectedSamplePoints, columnAvg, 0, cv::ReduceTypes::REDUCE_AVG); //dim=0 to reduce to single row

        projAvgPointsByCombo.push_back(columnAvg);

    }//end for each +/- combination

    ss << "The output, projAvgPointsByCombo: " << std::endl;
    ss << projAvgPointsByCombo << std::endl;

    //assemble 3 dimensional matrix of basis vectors with each of the testMultFactors applied
    //cv::Mat basisVectorsWithDifferentSigns(tempVectors.rows, tempVectors.cols, testMultFactors.rows, cv::DataType<double>::type);

    ////Do it with for loops until I am less sleepy

    //ss << "The size of basisVectorWithDifferentSigns is: " << tempVectors.rows << ", " << tempVectors.cols << ", " << testMultFactors.rows << "." << std::endl;

    //for (int combo = 0; combo < testMultFactors.rows; combo++) {
    //    for (int vec = 0; vec < tempVectors.cols; vec++) {
    //        double multFactor = (testMultFactors.at<int>(combo, vec) == 0) ? 1.0 : -1.0;
    //        for (int rgb = 0; rgb < tempVectors.rows; rgb++) {
    //            double tempVal = multFactor * tempVectors.at<double>(rgb, vec);

    //            ss << tempVal << std::endl;

    //            basisVectorsWithDifferentSigns.at<double>(vec, rgb, combo) = tempVal;
    //        }
    //    }
    //}

    //ss << "Did this work at all?" << std::endl;

    //for (int combo = 0; combo < testMultFactors.rows; combo++) {
    //    for (int vec = 0; vec < tempVectors.cols; vec++) {
    //        for (int rgb = 0; rgb < tempVectors.rows; rgb++) {
    //            ss << basisVectorsWithDifferentSigns.at<double>(vec, rgb, combo) << std::endl;
    //        }
    //    }
    //}

    //ss << basisVectorsWithDifferentSigns << std::endl;


    //Now do something with the subsampleOfPixels




    tempOut << ss.str() << std::endl;
    tempOut.close();



    //If the basis vectors were input as row vectors, transpose them back to that orientation

    //HACK! use the columnBasisVectors as the output, until I finish writing tests
    cv::Mat outputMatrix;
    if (basisVecDir == VectorDirection::COLUMNVECTORS) {
        //no change
        outputMatrix = columnBasisVectors.clone();
    }
    else if (basisVecDir == VectorDirection::ROWVECTORS) {
        //transpose the inputVectors matrix
        cv::transpose(columnBasisVectors, outputMatrix);
    }
    else {
        //invalid, but use the inputVectors as given and proceed anyway
        outputMatrix = columnBasisVectors.clone();
    }
    outputVectors.assign(outputMatrix);
}//end OptimizeBasisVectorSigns



void BasisTransform::CreatePixelSubsample(cv::InputArray sourcePixels, cv::OutputArray subsample, int numberOfPixels) {
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
