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

#ifndef SEDEEN_SRC_PLUGINS_STAINANALYSIS_ODCONVERSION_H
#define SEDEEN_SRC_PLUGINS_STAINANALYSIS_ODCONVERSION_H

#include <cmath>
#include <vector>

///A class with static methods or lookup-table conversion to/from optical density
class ODConversion {
public:
    ///Constructor to build the lookup table
    ODConversion() {
        //Build the lookup table
        m_convLookup.reserve(GetRGBMaxValue());
        for (int i = 0; i < GetRGBMaxValue(); i++) {
            m_convLookup.push_back(ConvertRGBtoOD(static_cast<double>(i)));
        }
    }//end lookup table constructor

    virtual ~ODConversion(void) {
        m_convLookup.clear();
        m_convLookup.shrink_to_fit();
    }//end destructor

    ///RGB to OD conversion using a lookup table
    inline const double LookupRGBtoOD(const int &_color) const {
        //Using 'at' instead of operator[] means an out_of_range exception can be thrown
        try {
            return m_convLookup.at(_color);
        }
        catch (const std::out_of_range&) {
            return ConvertRGBtoOD(static_cast<double>(_color));
        }
    }//end LookupRGBToOD

    ///OD to RGB conversion using a lookup table (imprecise - use carefully)
    inline const int LookupODtoRGB(const double &_OD) const {
        //Traverse the vector in reverse, find index at which _OD is smaller than stored value
        for (auto p = m_convLookup.rbegin(); p != m_convLookup.rend(); ++p) {
            if (_OD <= *p) {
                return static_cast<int>(m_convLookup.rend() - p);
            }
        }
        //if not found
        return static_cast<int>(ConvertODtoRGB(_OD));
    }//end LookupODToRGB

    ///Convert from color space (0 to 255 RGB value) to optical density
    inline static const double ConvertRGBtoOD(const double &_color) {
        double scaleMax = static_cast<double>(GetRGBMaxValue());
        //Avoid trying to calculate log(0)
        double color = (_color <= 0.0) ? GetODMinValue() : _color;
        double OD = (color == scaleMax) ? GetODMinValue() : -std::log10(color / scaleMax);
        //Push negative and 0 values up to small positive value
        OD = (OD < GetODMinValue()) ? GetODMinValue() : OD;
        return OD;
    }//end ConvertRGBtoOD

    ///Convert from optical density to color space (0 to 255 RGB value)
    inline static const double ConvertODtoRGB(const double &_OD) {
        double scaleMax = static_cast<double>(GetRGBMaxValue());
        //Push negative and 0 values up to small positive value
        double OD = (_OD < GetODMinValue()) ? GetODMinValue() : _OD;
        double color = std::round(scaleMax * std::pow(10.0, -OD));
        //Valid range is 0 to scaleMax
        color = (color < 0.0) ? 0.0 : color;
        color = (color > scaleMax) ? scaleMax : color;
        return color;
    }//end ConvertODtoRGB

public:
    ///Choose a value to represent near-zero in this class
    inline static const double GetODMinValue() { return 1e-6; }
    ///Define the maximum value of the RGB scale used in images
    inline static const int GetRGBMaxValue() { return 255; }

private:
    ///A lookup table implemented using a vector (relies on RGB values being integers)
    std::vector<double> m_convLookup;
    
};

#endif
