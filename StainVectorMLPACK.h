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
     
#ifndef SEDEEN_SRC_FILTER_STAINVECTORMLPACK_H
#define SEDEEN_SRC_FILTER_STAINVECTORMLPACK_H

#include "Global.h"
#include "Geometry.h"
#include "Image.h"

#include "StainVectorOpenCV.h"

//MLPACK includes
#include <mlpack/core.hpp>
#include <armadillo>

namespace sedeen {
namespace image {

class PATHCORE_IMAGE_API StainVectorMLPACK : public StainVectorOpenCV {
public:
    StainVectorMLPACK(std::shared_ptr<tile::Factory> source);
    ~StainVectorMLPACK();

    ///Utility method to check the equality of the contents of two Armadillo matrices (format used by MLPACK)
    template<class Ty>
    static const bool AreEqual(const arma::Mat<Ty> &array1, const arma::Mat<Ty> &array2) {
        // treat two empty arrays as identical
        if (array1.empty() && array2.empty()) {
            return true;
        }
        // if dimensionality is not identical, these arrays are not identical
        if (array1.n_cols != array2.n_cols || array1.n_rows != array2.n_rows || array1.size() != array2.size()) {
            return false;
        }
        //Use the Armadillo function approx_equal, which compares within a tolerance
        //Use reldiff to make this method more general, so that the hardcoded tolerance 
        //is relative to the scale of the largest element value
        Ty tolerance = static_cast<Ty>(1e-6); 
        return arma::approx_equal(array1, array2, "reldiff", tolerance);
    }//end AreEqual

    ///Convert an OpenCV matrix to an Armadillo matrix of templated type
    template<class Ty>
    static const arma::Mat<Ty> CVMatToArmaMat(cv::InputArray _input) {
        cv::Mat inputCVMat = _input.getMat();
        //Armadillo stores in column-major order, whereas OpenCV uses row-major, so transpose first
        cv::Mat transposedCV;
        cv::transpose(inputCVMat, transposedCV);
        arma::Mat<Ty> outArmaMat(reinterpret_cast<Ty*>(transposedCV.data), transposedCV.cols, transposedCV.rows);
        return outArmaMat;
    }//end CVMatToArmaMat

    ///Convert an Armadillo matrix to an OpenCV matrix (note: CV data is row-major order, Armadillo is column-major)
    template<class Ty>
    static cv::Mat ArmaMatToCVMat(const arma::Mat<Ty> &inputArmaMat) {
        int rows = static_cast<int>(inputArmaMat.n_rows);
        int cols = static_cast<int>(inputArmaMat.n_cols);
        cv::Mat typedMat = cv::Mat_<Ty>(cols, rows, const_cast<Ty*>(inputArmaMat.memptr()));
        cv::Mat transposedMat;// (typedMat);
        cv::transpose(typedMat, transposedMat);
        return transposedMat;
    }//end ArmaMatToCVMat

    ///Convert an Armadillo matrix and assign to a cv::OutputArray reference (note: CV data is row-major order, Armadillo is column-major)
    template<class Ty>
    static void ArmaMatToCVMat(const arma::Mat<Ty> &inputArmaMat, cv::OutputArray outArray) {
        outArray.assign(ArmaMatToCVMat<Ty>(inputArmaMat));
    }//end ArmaMatToCVMat

};

} // namespace image
} // namespace sedeen
#endif
