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

#ifndef STAINANALYSIS_BASISTRANSFORM_H
#define STAINANALYSIS_BASISTRANSFORM_H

#include <random>

//OpenCV include
#include <opencv2/core/core.hpp>

namespace sedeen {
namespace image {



class BasisTransform {
public:
    ///An enum to identify which axis input vectors are arranged in (outputs will match, if applicable)
    enum VectorDirection {
        COLUMNVECTORS,
        ROWVECTORS,
        UNDETERMINED
    };

public:
    BasisTransform();
    virtual ~BasisTransform();

    void ProjectionOnly(cv::InputArray sourcePoints, cv::OutputArray outputPoints, cv::InputArray basisVectors,
        const VectorDirection &sourcePointDir = VectorDirection::ROWVECTORS);



    ///Perform Principal Component Analysis (PCA) on a set of points, find a basis, transform points to new basis.
    void PCAPointTransform(cv::InputArray sourcePoints, cv::OutputArray outputPoints, 
        cv::InputArray sourceMask = cv::noArray(), cv::InputArray inputMean = cv::noArray(),
        const VectorDirection &sourcePointDir = VectorDirection::ROWVECTORS);

    ///Use the member variable basis vectors to create a set of points projected into a new basis. Set subtractMean to translate before projection.
    bool projectPoints(cv::InputArray sourcePoints, cv::OutputArray projectedPoints, bool subtractMean = true) const;

    ///Given a 2D projected point set, backproject to the original basis using the stored basis vectors. Set addMean to translate after back-projection.
    bool backProjectPoints(cv::InputArray projectedPoints, cv::OutputArray backProjPoints, bool addMean = true) const;

    ///Set/Get the numTestingPixels member variable
    inline void SetNumTestingPixels(const int n) { m_numTestingPixels = n; }
    ///Set/Get the numTestingPixels member variable
    inline const int GetNumTestingPixels() const { return m_numTestingPixels; }

    ///Get the basis vectors computed in this class, if not empty. Returns true on success, false if member Mat is empty.
    bool GetBasisVectors(cv::OutputArray basisVectors) const;
    ///Get the basis vectors computed in this class. Returns a possibly-empty matrix
    cv::Mat GetBasisVectors() const;

    ///Get the member point mean
    void GetPointMean(cv::OutputArray mean) const;
    ///Get the member point mean
    cv::Mat GetPointMean() const;
    ///Get some or all of the member eigenvalues (nVals = -1 to return all)
    void GetEigenvalues(cv::OutputArray evals, const int nVals = -1) const;
    ///Get some or all of the member eigenvalues (nVals = -1 to return all)
    cv::Mat GetEigenvalues(const int nVals = -1) const;
    ///Get the direction the eigenvectors are oriented in from the eigenvalues
    const VectorDirection GetEigenvectorElementsDirection() const;
    ///Get some or all of the member eigenvectors (nVecs = -1 to return all)
    void GetEigenvectors(cv::OutputArray evecs, const int nVecs = -1, 
        const VectorDirection &evecDir = VectorDirection::ROWVECTORS) const;
    ///Get some or all of the member eigenvectors (nVecs = -1 to return all)
    cv::Mat GetEigenvectors(const int nVecs = -1, 
        const VectorDirection &evecDir = VectorDirection::ROWVECTORS) const;




    ///Which signs should be used for the basis vectors? Test projecting some source points, try to get ++ quadrant projections
    void OptimizeBasisVectorSigns(cv::InputArray sourcePixels,
        cv::InputArray inputVectors, cv::OutputArray outputVectors, bool useMean = true,
        const VectorDirection &basisVecDir = VectorDirection::COLUMNVECTORS);



protected:
    ///Given points and a set of basis vectors, create a set of points projected into the new basis. Set subtractMean to translate before projection.
    void projectPoints(cv::InputArray sourcePoints, cv::OutputArray projectedPoints, 
        cv::InputArray basisVectors, cv::InputArray means, bool subtractMean = true) const;
    ///Given a 2D projected point set and a basis vectors, backproject to the original basis. Set addMean to translate after back-projection.
    void backProjectPoints(cv::InputArray projectedPoints, cv::OutputArray backProjPoints, 
        cv::InputArray basisVectors, cv::InputArray means, bool addMean = true) const;
    ///Randomly choose numberOfPixels rows from sourcePixels, copy them to the subsample OutputArray.
    void CreatePixelSubsample(cv::InputArray sourcePixels, cv::OutputArray subsample, const int numberOfPixels);
    ///Set the member basis vectors. Second parameter is direction of input vectors. Store basis vectors as row vectors.
    void SetBasisVectors(cv::InputArray basisVectors, const VectorDirection &vecDir = VectorDirection::ROWVECTORS);
    ///Set the member point mean
    void SetPointMean(cv::InputArray mean);
    ///Set the member eigenvalues
    void SetEigenvalues(cv::InputArray evals);
    ///Set the member eigenvectors
    void SetEigenvectors(cv::InputArray evecs);

protected:
    ///Allow derived classes access to the random number generator (64-bit Mersenne Twister)
    std::mt19937_64 m_rgen; //64-bit Mersenne Twister

private:
    //We specifically want two basis vectors
    const int m_reqdBasisVectors = 2;

    ///The number of pixels to test basis vector orientations with
    int m_numTestingPixels;

    ///Basis vectors
    cv::Mat m_basisVectors;
    ///The mean position to centre the data on
    cv::Mat m_pointMean;
    ///All eigenvalues of the covariance matrix, in descending order
    cv::Mat m_eigenvalues;
    ///All eigenvectors of the covariance matrix, in eigenvalue descending order
    cv::Mat m_eigenvectors;
};

} // namespace image
} // namespace sedeen
#endif
