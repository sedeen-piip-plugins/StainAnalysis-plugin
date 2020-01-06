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

#ifndef SEDEEN_SRC_PLUGINS_STAINANALYSIS_STAINVECTORMATH_H
#define SEDEEN_SRC_PLUGINS_STAINANALYSIS_STAINVECTORMATH_H

#include <vector>
#include <array>
#include <cmath>
#include <numeric>

///A class with static methods to operate on stain vectors
class StainVectorMath
{
public:
    ///Compute the inverse of a 3x3 matrix using Boost qvm
    static void Compute3x3MatrixInverse(const double (&inputMat)[9], double (&inversionMat)[9]);

    ///Make a 3x3 matrix, expressed as a 9-element array, have unitary rows. Preserve rows of zeros
    static void Make3x3MatrixUnitary(const double (&inputMat)[9], double (&unitaryMat)[9]);

    ///Check the Norm values of sets of three input elements, replace with a unitary row if zero
    static void ConvertZeroRowsToUnitary(const double (&inputMat)[9], double (&unitaryMat)[9]);

    ///Check the norm values of sets of elements, normalize the row of values in the third argument and replace zero rows with that
    static void ConvertZeroRowsToUnitary(const double (&inputMat)[9], double (&unitaryMat)[9], const double (&replacementVals)[3]);

    ///Check whether rows of the given matrix sum to zero, but do not have all zero values
    static std::array<bool, 3> RowSumZeroCheck(const double (&inputMat)[9]);

    ///Multiply a 3x3 matrix and a 3x1 vector to produce a 3x1 vector
    static void Multiply3x3MatrixAndVector(const double (&inputMat)[9], const double (&inputVec)[3], double (&outputVec)[3]);

    ///Return an array of values of type Ty with size N normalized to unit length. Returns input array if norm is 0.
    template<class Ty, std::size_t N> 
    static std::array<Ty, N> NormalizeArray(std::array<Ty, N> arr) {
        std::array<Ty, N> out;
        Ty norm = Norm<std::array<Ty, N>::iterator, Ty>(arr.begin(), arr.end());
        //Check if the norm is zero. Return the input array if so.
        //Compare against C++11 zero initialization of the type Ty
        //Also check if the input container is empty
        if ((norm == Ty{}) || (arr.empty())) {
            return arr;
        }
        else {
            //Copy the input array to the out array
            std::copy(arr.begin(), arr.end(), out.begin());
            //Iterate through the out array, divide values by norm
            for (auto p = out.begin(); p != out.end(); ++p) {
                *p = static_cast<Ty>(*p / norm);
            }
            return out;
        }
    }//end NormalizeArray

    ///Calculate the norm of all the elements in a container, where each element is of type Ty
    template<typename Iter_T, class Ty>
    static Ty Norm(Iter_T first, Iter_T last) {
        return static_cast<Ty>(sqrt(std::inner_product(first, last, first, Ty{}))); //Use C++11 zero initialization of type Ty
    }//end Norm
};

#endif
