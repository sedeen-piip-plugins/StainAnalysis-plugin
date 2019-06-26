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

#include "StainVectorMath.h"

//Boost includes
#include <boost/qvm/vec.hpp>
#include <boost/qvm/mat.hpp>
#include <boost/qvm/vec_access.hpp>
#include <boost/qvm/mat_access.hpp>
#include <boost/qvm/vec_operations.hpp>
#include <boost/qvm/mat_operations.hpp>
#include <boost/qvm/vec_mat_operations.hpp>
#include <boost/qvm/map_mat_mat.hpp>

///Compute the inverse of a 3x3 matrix using Boost qvm: ensure matrix is unitary before using
void StainVectorMath::Compute3x3MatrixInverse(double inputMat[9], double inversionMat[9]) {
    //Clear the inversionMat (output) array
    for (int i = 0; i < 9; i++) { inversionMat[i] = 0.0; }
    //Define an output matrix
    boost::qvm::mat<double, 3, 3> outputMatrix;

    //Reshape to a 3x3 matrix
    double reshapedInput[3][3] = { 0.0 };
    int i, j = 0;
    for (int x = 0; x < 9; x++) {
        i = x / 3; //integer division
        j = x % 3;
        reshapedInput[i][j] = inputMat[x];
    }
     
    //Check the determinant of the matrix: is it 0? Thus is the matrix singular?
    double theDeterminant = boost::qvm::determinant(reshapedInput);
    //Get the inverse of the reshapedInput matrix if the determinant is not zero
    //Return matrix of all zeros if the determinant is zero
    if (abs(theDeterminant) < StainVectorMath::GetODMinValue()) {
        outputMatrix = boost::qvm::zero_mat<double, 3, 3>();
    }
    else {
        outputMatrix = boost::qvm::transposed(boost::qvm::inverse(reshapedInput));
    }

    //Assign output values (first index is row), no looping
    inversionMat[0] = boost::qvm::A<0, 0>(outputMatrix);
    inversionMat[1] = boost::qvm::A<0, 1>(outputMatrix);
    inversionMat[2] = boost::qvm::A<0, 2>(outputMatrix);
    inversionMat[3] = boost::qvm::A<1, 0>(outputMatrix);
    inversionMat[4] = boost::qvm::A<1, 1>(outputMatrix);
    inversionMat[5] = boost::qvm::A<1, 2>(outputMatrix);
    inversionMat[6] = boost::qvm::A<2, 0>(outputMatrix);
    inversionMat[7] = boost::qvm::A<2, 1>(outputMatrix);
    inversionMat[8] = boost::qvm::A<2, 2>(outputMatrix);
    //void return
}//end Compute3x3MatrixInverse

void StainVectorMath::Make3x3MatrixUnitary(double inputMat[9], double unitaryMat[9]) {
    //Bundle the input values in rows of three
    std::vector<std::array<double, 3>> inputRows;
    inputRows.push_back(std::array<double, 3>({ inputMat[0], inputMat[1], inputMat[2] }));
    inputRows.push_back(std::array<double, 3>({ inputMat[3], inputMat[4], inputMat[5] }));
    inputRows.push_back(std::array<double, 3>({ inputMat[6], inputMat[7], inputMat[8] }));
    //Get the norm values 
    std::vector<double> normVals;
    normVals.push_back(StainVectorMath::Norm<std::array<double, 3>::iterator, double>(inputRows[0].begin(), inputRows[0].end()));
    normVals.push_back(StainVectorMath::Norm<std::array<double, 3>::iterator, double>(inputRows[1].begin(), inputRows[1].end()));
    normVals.push_back(StainVectorMath::Norm<std::array<double, 3>::iterator, double>(inputRows[2].begin(), inputRows[2].end()));
    //Create output rows
    std::vector<std::array<double, 3>> outputRows;
    for (auto it = normVals.begin(); it != normVals.end(); ++it) {
        bool smallValue = (*it < 10.0*StainVectorMath::GetODMinValue());
        std::array<double, 3> addRow = inputRows[it - normVals.begin()];
        if (smallValue) {
            //do not modify
        }
        else {
            for (auto p = addRow.begin(); p != addRow.end(); ++p) {
                *p = *p / *it;
            }
        }
        outputRows.push_back(addRow);
    }
    //Assign to the unitary (output) matrix
    for (int x = 0; x < 9; x++) {
        int i = x / 3;
        int j = x % 3;
        unitaryMat[x] = outputRows[i][j];
    }
}//end Make3x3MatrixUnitary

void StainVectorMath::ConvertZeroRowsToUnitary(double inputMat[9], double unitaryMat[9]) {
    double replacementVals[3] = { 1.0,1.0,1.0 };
    StainVectorMath::ConvertZeroRowsToUnitary(inputMat, unitaryMat, replacementVals);
}//end ConvertZeroRowsToUnitary

void StainVectorMath::ConvertZeroRowsToUnitary(double inputMat[9], double unitaryMat[9], double replacementVals[3]) {
    //Bundle the input values in rows of three
    std::vector<std::array<double, 3>> inputRows;
    inputRows.push_back(std::array<double, 3>({ inputMat[0], inputMat[1], inputMat[2] }));
    inputRows.push_back(std::array<double, 3>({ inputMat[3], inputMat[4], inputMat[5] }));
    inputRows.push_back(std::array<double, 3>({ inputMat[6], inputMat[7], inputMat[8] }));
    //Get the norm values 
    std::vector<double> normVals;
    normVals.push_back(StainVectorMath::Norm<std::array<double, 3>::iterator, double>(inputRows[0].begin(), inputRows[0].end()));
    normVals.push_back(StainVectorMath::Norm<std::array<double, 3>::iterator, double>(inputRows[1].begin(), inputRows[1].end()));
    normVals.push_back(StainVectorMath::Norm<std::array<double, 3>::iterator, double>(inputRows[2].begin(), inputRows[2].end()));
    //Compare against 10xGetODMinValue, set to unitary row if smaller
    std::array<double, 3> replacementArray = { replacementVals[0], replacementVals[1], replacementVals[2] };
    auto unitaryRow = StainVectorMath::NormalizeArray(replacementArray);
    std::vector<std::array<double, 3>> outputRows;
    for (auto it = normVals.begin(); it != normVals.end(); ++it) {
        bool smallValue = (*it < 10.0*StainVectorMath::GetODMinValue());
        auto addRow = smallValue ? unitaryRow : inputRows[it-normVals.begin()];
        outputRows.push_back(addRow);
    }
    //Assign to the unitary (output) matrix
    for (int x = 0; x < 9; x++) {
        int i = x / 3;
        int j = x % 3;
        unitaryMat[x] = outputRows[i][j];
    }
}//end ConvertZeroRowsToUnitary

std::array<bool, 3> StainVectorMath::RowSumZeroCheck(double inputMat[9]) {
    std::array<bool, 3> returnVals;
    //Bundle the input values in rows of three
    std::vector<std::array<double, 3>> inputRows;
    inputRows.push_back(std::array<double, 3>({ inputMat[0], inputMat[1], inputMat[2] }));
    inputRows.push_back(std::array<double, 3>({ inputMat[3], inputMat[4], inputMat[5] }));
    inputRows.push_back(std::array<double, 3>({ inputMat[6], inputMat[7], inputMat[8] }));
    //Get the sums
    std::vector<double> rowSums;
    rowSums.push_back(std::accumulate(inputRows[0].begin(), inputRows[0].end(), 0.0));
    rowSums.push_back(std::accumulate(inputRows[1].begin(), inputRows[1].end(), 0.0));
    rowSums.push_back(std::accumulate(inputRows[2].begin(), inputRows[2].end(), 0.0));
    //Get the norm values 
    std::vector<double> normVals;
    normVals.push_back(StainVectorMath::Norm<std::array<double, 3>::iterator, double>(inputRows[0].begin(), inputRows[0].end()));
    normVals.push_back(StainVectorMath::Norm<std::array<double, 3>::iterator, double>(inputRows[1].begin(), inputRows[1].end()));
    normVals.push_back(StainVectorMath::Norm<std::array<double, 3>::iterator, double>(inputRows[2].begin(), inputRows[2].end()));

    //For each, if rowSums is zero and normVals is non-zero, returnVals is true
    for (int i = 0; i < 3; i++) {
        bool checkVal = ((abs(rowSums[i]) < StainVectorMath::GetODMinValue()) && (normVals[i] > 0.0));
        returnVals[i] = checkVal ? true : false;
    }
    return returnVals;
}//end RowSumZeroCheck

void StainVectorMath::Multiply3x3MatrixAndVector(double inputMat[9], double inputVec[3], double outputVec[3]) {
    //Clear the output vector
    for (int i = 0; i < 3; i++) { outputVec[i] = 0.0; }
    //Reshape to a 3x3 matrix
    double reshapedMatrix[3][3] = { 0.0 };
    int i, j = 0;
    for (int x = 0; x < 9; x++) {
        i = x / 3;
        j = x % 3;
        reshapedMatrix[i][j] = inputMat[x];
    }

    //Assign matrix and input vector to boost::qvm types
    boost::qvm::mat<double, 3, 3> inputQVMMatrix = boost::qvm::mref(reshapedMatrix);
    boost::qvm::vec<double, 3> inputQVMVector = boost::qvm::zero_vec<double, 3>();
    boost::qvm::A<0>(inputQVMVector) = inputVec[0];
    boost::qvm::A<1>(inputQVMVector) = inputVec[1];
    boost::qvm::A<2>(inputQVMVector) = inputVec[2];

    //Multiplication
    boost::qvm::vec<double, 3> outputQVMVector;
    outputQVMVector = inputQVMMatrix * inputQVMVector;

    //Assign to the output vector
    outputVec[0] = boost::qvm::A<0>(outputQVMVector);
    outputVec[1] = boost::qvm::A<1>(outputQVMVector);
    outputVec[2] = boost::qvm::A<2>(outputQVMVector);
    //void return
}//end Multiply3x3MatrixAndVector
