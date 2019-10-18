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

#ifndef STAINANALYSIS_BASISTRANSFORM_H
#define STAINANALYSIS_BASISTRANSFORM_H

#include <cassert>
#include <random>

//OpenCV include
#include <opencv2/core/core.hpp>
//#include <opencv2/imgproc/imgproc.hpp>
//#include <opencv2/highgui/highgui.hpp>


namespace sedeen {
namespace image {

class BasisTransform {
public:
    ///An enum to identify which axis the input vectors are arranged in (outputs will match)
    enum VectorDirection {
        COLUMNVECTORS,
        ROWVECTORS
    };

public:
    BasisTransform();
    ~BasisTransform();

    ///Perform Principal Component Analysis (PCA) on a set of points, find a basis, transform points to new basis.
    void PCAPointTransform(cv::InputArray sourcePoints, cv::OutputArray outputPoints);

    ///Which signs should be used for the basis vectors? Test projecting some source points, try to get ++ quadrant projections
    void OptimizeBasisVectorSigns(cv::InputArray sourcePixels, 
        cv::InputArray inputVectors, cv::OutputArray outputVectors, 
        VectorDirection basisVecDir = VectorDirection::COLUMNVECTORS);

    ///Given points and a set of basis vectors, create a set of points projected into the new basis
    void projectPoints(cv::InputArray sourcePoints, cv::InputArray basisVectors, cv::OutputArray projectedPoints);

    ///Given points that have been projected into a given basis, backproject to the original basis
    void backProjectPoints(cv::InputArray projectedPoints, cv::InputArray basisVectors, cv::OutputArray backProjPoints);

    ///Get the basis vectors computed in this class, or an empty matrix if they have not been computed yet.
    void GetBasisVectors(cv::OutputArray basisVectors);

    ///Set/Get the numTestingPixels member variable
    inline void SetNumTestingPixels(int n) { m_numTestingPixels = n; }
    ///Set/Get the numTestingPixels member variable
    inline  int GetNumTestingPixels() { return m_numTestingPixels; }

protected:
    ///Randomly choose numberOfPixels rows from sourcePixels, copy them to the subsample OutputArray.
    void CreatePixelSubsample(cv::InputArray sourcePixels, cv::OutputArray subsample, int numberOfPixels);


protected:
    ///Allow derived classes access to the random number generator (64-bit Mersenne Twister)
    std::mt19937_64 m_rgen; //64-bit Mersenne Twister

private:
    int m_numTestingPixels;

};

} // namespace image
} // namespace sedeen
#endif
