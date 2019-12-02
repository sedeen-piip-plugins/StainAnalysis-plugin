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

#include "StainVectorOpenCV.h"
#include "StainVectorMath.h"

namespace sedeen {
namespace image {

StainVectorOpenCV::StainVectorOpenCV(std::shared_ptr<tile::Factory> source) 
    : StainVectorBase(source)
{}//end constructor

StainVectorOpenCV::~StainVectorOpenCV(void) {
}//end destructor

void StainVectorOpenCV::StainCVMatToCArray(cv::InputArray inputData, double (&outputVectors)[9], bool normalize /*=false*/) {
    double tempOutput[9] = { 0.0 };
    if (inputData.empty()) { return; }
    cv::Mat inputMatSquare, _inputMat(inputData.getMat());
    _inputMat.convertTo(inputMatSquare, cv::DataType<double>::type);
    //Reshape the matrix and get the data
    int numElements = static_cast<int>(inputMatSquare.total());
    if (numElements % 3 != 0) { return; }
    //Reshape the matrix (does not reallocate memory)
    cv::Mat inputMatFlat = inputMatSquare.reshape(0, 1);
    //Fill the tempOutput array with data elements
    std::copy(inputMatFlat.begin<double>(), inputMatFlat.end<double>(), std::begin(tempOutput));
    //for (auto p = inputMatFlat.begin<double>(); p != inputMatFlat.end<double>(); p++) {
    //    int el = static_cast<int>(p.lpos());
    //    tempOutput[el] = *p;
    //}
    //Pad the tempOutput array if necessary
    std::fill(std::begin(tempOutput) + numElements, std::end(tempOutput), 0.0);
    //for (int q = numElements; q < 9; q++) {
    //    tempOutput[q] = 0.0;
    //}
    //If normalize is true, make outputVectors normalized, else directly copy from tempOutput
    if (normalize) {
        StainVectorMath::Make3x3MatrixUnitary(tempOutput, outputVectors);
    }
    else {
        std::copy(std::begin(tempOutput), std::end(tempOutput), std::begin(outputVectors));
        //for (int i = 0; i < 9; i++) {
        //    outputVectors[i] = tempOutput[i];
        //}
    }
}//end StainCVMatToCArray

void StainVectorOpenCV::StainCArrayToCVMat(double (&inputVectors)[9], cv::OutputArray outputData, bool normalize /*=false*/) {
    //If normalize is true, fill an array made unitary, else make a copy of the input data
    double inputCopy[9] = { 0.0 };
    if (normalize) {
        StainVectorMath::Make3x3MatrixUnitary(inputVectors, inputCopy);
    }
    else {
        std::copy(std::begin(inputVectors), std::end(inputVectors), std::begin(inputCopy));
        //for (int i = 0; i < 9; i++) {
        //    inputCopy[i] = inputVectors[i];
        //}
    }
    //Create a cv::Mat of type double, copy inputCopy data
    cv::Mat inputMatFlat(1, 9, cv::DataType<double>::type);
    std::copy(std::begin(inputCopy), std::end(inputCopy), inputMatFlat.begin<double>());
    //for (auto p = inputMatFlat.begin<double>(); p != inputMatFlat.end<double>(); p++) {
    //    int el = static_cast<int>(p.lpos());
    //    *p = inputCopy[el];
    //}
    //Reshape the matrix (does not reallocate or copy data)
    cv::Mat inputMatSquare = inputMatFlat.reshape(0, 3);
    outputData.assign(inputMatSquare);
}//end StainCArrayToCVMat

} // namespace image
} // namespace sedeen
