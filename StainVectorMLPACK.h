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
//#include <armadillo>


namespace sedeen {
namespace image {

class PATHCORE_IMAGE_API StainVectorMLPACK : public StainVectorOpenCV {
public:
    StainVectorMLPACK(std::shared_ptr<tile::Factory> source);
    ~StainVectorMLPACK();

    ///Utility method to check the equality of the contents of two Armadillo matrices (format used by MLPACK)
    bool AreEqual(cv::InputArray array1, cv::InputArray array2);

protected:
    ///Convert stain vector data as 9-element C array to Armadillo matrix (as row vectors)
    void StainCArrayToArmaMat(double (&inutVectors)[9], cv::OutputArray outputData, 
        const bool normalize = false, const int _numRows = -1);
    ///Convert stain vector data from Armadillo matrix (as row vectors) to 9-element C array
    void StainArmaMatToCArray(cv::InputArray inputData, double (&outputVectors)[9], const bool normalize = false);

    ///Convert an OpenCV matrix to an Armadillo matrix (note: CV is row-major order, Armadillo is column-major)



    ///Convert an Armadillo matrix to an OpenCV matrix (note: CV is row-major order, Armadillo is column-major)



};

} // namespace image
} // namespace sedeen
#endif
