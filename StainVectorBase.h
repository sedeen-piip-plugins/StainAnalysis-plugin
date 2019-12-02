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

#ifndef SEDEEN_SRC_FILTER_STAINVECTORBASE_H
#define SEDEEN_SRC_FILTER_STAINVECTORBASE_H

#include "Global.h"
#include "Geometry.h"
#include "Image.h"

#include "RandomWSISampler.h"

namespace sedeen {
namespace image {

class PATHCORE_IMAGE_API StainVectorBase {
public:
    StainVectorBase(std::shared_ptr<tile::Factory> source);
    ~StainVectorBase();

    ///The core functionality of a stain vector class; fills the 9-element array with three stain vectors
    virtual void ComputeStainVectors(double (&outputVectors)[9]);

protected:
    ///Returns a shared pointer to the source factory, protected so only derived classes may access it
    inline std::shared_ptr<tile::Factory> GetSourceFactory() { return m_sourceFactory; }
    ///Sets the source factory, protected so only derived classes may modify it
    inline void SetSourceFactory(std::shared_ptr<tile::Factory> source) { m_sourceFactory = source; }
    ///Access the random pixel chooser
    inline std::shared_ptr<RandomWSISampler> GetRandomWSISampler() { return m_randomWSISampler; }

private:
    std::shared_ptr<tile::Factory> m_sourceFactory;
    std::shared_ptr<RandomWSISampler> m_randomWSISampler;
};

} // namespace image
} // namespace sedeen
#endif
