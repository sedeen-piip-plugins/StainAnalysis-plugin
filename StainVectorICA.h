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

#ifndef SEDEEN_SRC_FILTER_STAINVECTORICA_H
#define SEDEEN_SRC_FILTER_STAINVECTORICA_H

#include "Global.h"
#include "Geometry.h"
#include "Image.h"

#include "StainVectorMLPACK.h"

namespace sedeen {
namespace image {

class PATHCORE_IMAGE_API StainVectorICA : public StainVectorMLPACK {
public:
    StainVectorICA(std::shared_ptr<tile::Factory> source, double ODthreshold = 0.15);
    ~StainVectorICA();

    ///Fill the 9-element array with three stain vectors
    virtual void ComputeStainVectors(double(&outputVectors)[9]);
    ///Overload of the basic method, includes sampleSize parameter
    void ComputeStainVectors(double(&outputVectors)[9], const int sampleSize);

    ///Get/Set the average optical density threshold
    inline const double GetODThreshold() const { return m_avgODThreshold; }
    ///Get/Set the average optical density threshold
    inline void SetODThreshold(const double t) { m_avgODThreshold = t; }

    ///Get/Set the sample size, the number of pixels to choose
    inline const int GetSampleSize() const { return m_sampleSize; }
    ///Get/Set the sample size, the number of pixels to choose
    inline void SetSampleSize(const int s) { m_sampleSize = s; }

protected:
    ///Get/Set the number of stains
    inline const int GetNumStains() const { return m_numStains; }
    ///Get/Set the number of stains
    inline void SetNumStains(const int n) { m_numStains = n; }

private:
    double m_avgODThreshold;

    ///The number of pixels that should be used to calculate the stain vectors
    int m_sampleSize;
    ///The number of stains to obtain. Set to be 2 in member initialization of constructor.
    int m_numStains;
};

} // namespace image
} // namespace sedeen
#endif
