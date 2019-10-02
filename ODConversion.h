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

#ifndef SEDEEN_SRC_PLUGINS_STAINANALYSIS_ODCONVERSION_H
#define SEDEEN_SRC_PLUGINS_STAINANALYSIS_ODCONVERSION_H

#include <cmath>

///A class with static methods to operate on stain vectors
class ODConversion
{
public:
    ///Convert from color space (0 to 255 RGB value) to optical density
    inline static const double ConvertRGBtoOD(double color) {
        double scaleMax = 255.0;
        //Avoid trying to calculate log(0)
        color = (color <= 0.0) ? ODConversion::GetODMinValue() : color;
        double OD = -std::log10(color / scaleMax);
        //Push negative and 0 values up to small positive value
        OD = (OD < ODConversion::GetODMinValue()) ? ODConversion::GetODMinValue() : OD;
        return OD;
    }//end ConvertRGBtoOD

    ///Convert from optical density to color space (0 to 255 RGB value)
    inline static const double ConvertODtoRGB(double OD) {
        double scaleMax = 255.0;
        //Push negative and 0 values up to small positive value
        OD = (OD < ODConversion::GetODMinValue()) ? ODConversion::GetODMinValue() : OD;
        double color = std::round(scaleMax * std::pow(10.0, -OD)); //rounds values away from 0
        //Valid range is 0 to scaleMax
        color = (color < 0.0) ? 0.0 : color;
        color = (color > 255.0) ? 255.0 : color;
        return color;
    }//end ConvertODtoRGB

public:
    ///Choose a value to represent near-zero in this class
    inline static const double GetODMinValue() { return 1e-6; }
};

#endif