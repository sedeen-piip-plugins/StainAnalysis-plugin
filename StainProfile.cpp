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

#include "StainProfile.h"

#include <stdexcept>
#include <cmath>
#include <sstream>
#include <fstream>
#include <iterator>
#include <tinyxml2.h>

StainProfile::StainProfile() 
    : m_stainAnalysisModelOptions( { "Ruifrok+Johnston Deconvolution" } ),
    m_stainSeparationAlgorithmOptions( {"Region-of-Interest Selection", 
        "Macenko Decomposition", "Non-Negative Matrix Factorization" } )
{
    //Build the XML document structure
    BuildXMLDocument();
}//end constructor

StainProfile::StainProfile(StainProfile &s) 
    : StainProfile()  //invoking the no-parameter constructor will build the XML document
{
    this->SetNameOfStainProfile(s.GetNameOfStainProfile());
    this->SetNumberOfStainComponents(s.GetNumberOfStainComponents());
    this->SetNameOfStainOne(s.GetNameOfStainOne());
    this->SetNameOfStainTwo(s.GetNameOfStainTwo());
    this->SetNameOfStainThree(s.GetNameOfStainThree());
    this->SetNameOfStainAnalysisModel(s.GetNameOfStainAnalysisModel());
    this->SetNameOfStainSeparationAlgorithm(s.GetNameOfStainSeparationAlgorithm());
    this->SetAllAnalysisModelParameters(s.GetAllAnalysisModelParameters());
    this->SetAllSeparationAlgorithmParameters(s.GetAllSeparationAlgorithmParameters());
    this->SetStainOneRGB(s.GetStainOneRGB());
    this->SetStainTwoRGB(s.GetStainTwoRGB());
    this->SetStainThreeRGB(s.GetStainThreeRGB());
}//end copy constructor

StainProfile::~StainProfile() {
    //Smartpointer used for m_xmlDoc, 
    //and an XMLDocument handles memory management for all its child objects.
    //No explicit delete required.
}//end destructor

bool StainProfile::SetNameOfStainProfile(const std::string &name) {
    //Direct assignment. Add checks if necessary.
    if (m_rootElement == nullptr) {
        return false;
    }
    else {
        m_rootElement->SetAttribute(nameOfStainProfileAttribute(), name.c_str());
        return true;
    }
}//end SetNameOfStainProfile

const std::string StainProfile::GetNameOfStainProfile() const {
    std::string result;
    if (m_rootElement == nullptr) {
        result = std::string();
    }
    else {
        const char* name = m_rootElement->Attribute(nameOfStainProfileAttribute());
        result = (name == nullptr) ? std::string() : std::string(name);
    }
    return result;
}//end GetNameOfStainProfile

bool StainProfile::SetNumberOfStainComponents(const int &components) {
    if (m_componentsElement == nullptr) {
        return false;
    }
    //Check if number is 0 or greater.
    //If not, do not assign value and return false
    if (components >= 0) {
        m_componentsElement->SetAttribute(numberOfStainsAttribute(), components);
        return true;
    }
    else {
        m_componentsElement->SetAttribute(numberOfStainsAttribute(), -1);
        return false;
    }
}//end SetNumberOfStainComponents

const int StainProfile::GetNumberOfStainComponents() const {
    int components(-1);
    if (m_componentsElement == nullptr) {
        return -1;
    }
    else {
        tinyxml2::XMLError eResult = m_componentsElement->QueryIntAttribute(numberOfStainsAttribute(), &components);
        if (eResult != tinyxml2::XML_SUCCESS) {
            return -1;
        }
    }
    return components;
}//end GetNumberOfStainComponents

bool StainProfile::SetNameOfStainOne(const std::string &name) {
    //Direct assignment.
    if (m_stainOneElement == nullptr) {
        return false;
    }
    else {
        m_stainOneElement->SetAttribute(nameOfStainAttribute(), name.c_str());
        return true;
    }
}//end SetNameOfStainOne

const std::string StainProfile::GetNameOfStainOne() const {
    std::string result;
    if (m_stainOneElement == nullptr) {
        result = std::string();
    }
    else {
        const char* name = m_stainOneElement->Attribute(nameOfStainAttribute());
        result = (name == nullptr) ? std::string() : std::string(name);
    }
    return result;
}//end GetNameOfStainOne

bool StainProfile::SetNameOfStainTwo(const std::string &name) {
    //Direct assignment.
    if (m_stainTwoElement == nullptr) {
        return false;
    }
    else {
        m_stainTwoElement->SetAttribute(nameOfStainAttribute(), name.c_str());
        return true;
    }
}//end SetNameOfStainTwo

const std::string StainProfile::GetNameOfStainTwo() const {
    std::string result;
    if (m_stainTwoElement == nullptr) {
        result = std::string();
    }
    else {
        const char* name = m_stainTwoElement->Attribute(nameOfStainAttribute());
        result = (name == nullptr) ? std::string() : std::string(name);
    }
    return result;
}//end GetNameOfStainTwo

bool StainProfile::SetNameOfStainThree(const std::string &name) {
    //Direct assignment
    if (m_stainThreeElement == nullptr) {
        return false;
    }
    else {
        m_stainThreeElement->SetAttribute(nameOfStainAttribute(), name.c_str());
        return true;
    }
}//end SetNameOfStainThree

const std::string StainProfile::GetNameOfStainThree() const {
    std::string result;
    if (m_stainThreeElement == nullptr) {
        result = std::string();
    }
    else {
        const char* name = m_stainThreeElement->Attribute(nameOfStainAttribute());
        result = (name == nullptr) ? std::string() : std::string(name);
    }
    return result;
}//end GetNameOfStainThree

bool StainProfile::SetNameOfStainAnalysisModel(const std::string &name) {
    //Check if the stain analysis model name given is in the valid list
    //If not, do not assign value and return false
    if (std::find(m_stainAnalysisModelOptions.begin(), m_stainAnalysisModelOptions.end(), name)
        != m_stainAnalysisModelOptions.end()) {
        //Found. Assign name 
        m_analysisModelElement->SetAttribute(analysisModelNameAttribute(), name.c_str());
        return true;
    }//else
    return false;
}//end SetNameOfStainAnalysisModel

const std::string StainProfile::GetNameOfStainAnalysisModel() const {
    std::string result;
    if (m_analysisModelElement == nullptr) {
        result = std::string();
    }
    else {
        const char* am = m_analysisModelElement->Attribute(analysisModelNameAttribute());
        result = (am == nullptr) ? std::string() : std::string(am);
    }
    return result;
}//end GetNameOfStainAnalysisModel

bool StainProfile::SetNameOfStainSeparationAlgorithm(const std::string &name) {
    //Check if the stain separation algorithm name given is in the valid list
    //If not, do not assign value and return false
    if (std::find(m_stainSeparationAlgorithmOptions.begin(), m_stainSeparationAlgorithmOptions.end(), name) 
        != m_stainSeparationAlgorithmOptions.end()) {
        //Found. Assign name 
        m_algorithmElement->SetAttribute(algorithmNameAttribute(), name.c_str());
        return true;
    }//else
    return false;
}//end SetNameOfStainSeparationAlgorithm

const std::string StainProfile::GetNameOfStainSeparationAlgorithm() const {
    std::string result;
    if (m_algorithmElement == nullptr) {
        result = std::string();
    }
    else {
        const char* alg = m_algorithmElement->Attribute(algorithmNameAttribute());
        result = (alg == nullptr) ? std::string() : std::string(alg);
    }
    return result;
}//end GetNameOfStainSeparationAlgorithm

bool StainProfile::SetStainOneRGB(const double &r, const double &g, const double &b) {
    //Place data into a std::array, pass to other overload of this method
    std::array<double, 3> rgb = { r,g,b };
    return SetStainOneRGB(rgb);
}//end SetStainOneRGB

bool StainProfile::SetStainOneRGB(double rgb[])
{
    //Check length of rgb first
    //This is a C-style array, so get the size the old way
    int len = sizeof(rgb) / sizeof(*rgb); //(total size divided by size of first element)
    if (len != 3) {
        return false;
    }
    else {
        //Place data into a std::array, pass to other overload of this method
        std::array<double, 3> arrayRGB;
        for (int i = 0; i < 3; ++i) {
            try {
                arrayRGB.at(i) = rgb[i];
            }
            catch (const std::out_of_range& rangeerr) {
                rangeerr.what();
                //The index is out of range. Return false
                return false;
            }
        }
        return SetStainOneRGB(arrayRGB);
    }
    return false;
}//end SetStainOneRGB(C array)

bool StainProfile::SetStainOneRGB(const std::array<double, 3> &rgb_in) {
    //Normalize the array before assigning
    std::array<double, 3> rgb = StainVectorMath::NormalizeArray<double, 3>(rgb_in);
    //Get the first stain value element, or return false if not found
    if (m_stainOneElement == nullptr) { return false; }
    //else
    tinyxml2::XMLElement* sVals = m_stainOneElement->FirstChildElement(stainValueTag());
    if (sVals == nullptr) { return false; }
    while (sVals != nullptr) {
        if (sVals->Attribute(valueTypeAttribute(), "r")) {
            sVals->SetText(rgb[0]);
        }
        else if (sVals->Attribute(valueTypeAttribute(), "g")) {
            sVals->SetText(rgb[1]);
        }
        else if (sVals->Attribute(valueTypeAttribute(), "b")) {
            sVals->SetText(rgb[2]);
        }
        sVals = sVals->NextSiblingElement(stainValueTag());
    }
    return true;
}//end SetStainOneRGB(std::array)

const std::array<double, 3> StainProfile::GetStainOneRGB() const {
    //Create a zero array to return
    std::array<double, 3> out = { 0.0,0.0,0.0 };
    double r, g, b; //temp values
    //Get the first stain value element, or return false if not found
    if (m_stainOneElement == nullptr) { return out; }
    tinyxml2::XMLElement* sVals = m_stainOneElement->FirstChildElement(stainValueTag());
    if (sVals == nullptr) { return out; }
    while (sVals != nullptr) {
        if (sVals->Attribute(valueTypeAttribute(), "r")) {
            tinyxml2::XMLError rResult = sVals->QueryDoubleText(&r);
            if (rResult != tinyxml2::XML_SUCCESS) { return out; }
        }
        else if (sVals->Attribute(valueTypeAttribute(), "g")) {
            tinyxml2::XMLError gResult = sVals->QueryDoubleText(&g);
            if (gResult != tinyxml2::XML_SUCCESS) { return out; }
        }
        else if (sVals->Attribute(valueTypeAttribute(), "b")) {
            tinyxml2::XMLError bResult = sVals->QueryDoubleText(&b);
            if (bResult != tinyxml2::XML_SUCCESS) { return out; }
        }
        sVals = sVals->NextSiblingElement(stainValueTag());
    }
    out[0] = r;
    out[1] = g;
    out[2] = b;
    return out;
}//end GetStainOneRGB

bool StainProfile::SetStainTwoRGB(const double &r, const double &g, const double &b) {
    //Place data into a std::array, pass to other overload of this method
    std::array<double, 3> rgb = { r,g,b };
    return SetStainTwoRGB(rgb);
}//end SetStainTwoRGB

bool StainProfile::SetStainTwoRGB(double rgb[]) {
    //Check length of rgb first
    //This is a C-style array, so get the size the old way
    int len = sizeof(rgb) / sizeof(*rgb);
    if (len != 3) {
        return false;
    }
    else {
        //Place data into a std::array, pass to other overload of this method
        std::array<double, 3> arrayRGB;
        for (int i = 0; i < 3; ++i) {
            try {
                arrayRGB.at(i) = rgb[i];
            }
            catch (const std::out_of_range& rangeerr) {
                rangeerr.what();
                //The index is out of range. Return false
                return false;
            }
        }
        return SetStainTwoRGB(arrayRGB);
    }
    return false;
}//end SetStainTwoRGB(C array)

bool StainProfile::SetStainTwoRGB(const std::array<double, 3> &rgb_in) {
    //Normalize the array before assigning
    std::array<double, 3> rgb = StainVectorMath::NormalizeArray<double, 3>(rgb_in);
    //Get the first stain value element, or return false if not found
    if (m_stainTwoElement == nullptr) { return false; }
    tinyxml2::XMLElement* sVals = m_stainTwoElement->FirstChildElement(stainValueTag());
    if (sVals == nullptr) { return false; }
    while (sVals != nullptr) {
        if (sVals->Attribute(valueTypeAttribute(), "r")) {
            sVals->SetText(rgb[0]);
        }
        else if (sVals->Attribute(valueTypeAttribute(), "g")) {
            sVals->SetText(rgb[1]);
        }
        else if (sVals->Attribute(valueTypeAttribute(), "b")) {
            sVals->SetText(rgb[2]);
        }
        sVals = sVals->NextSiblingElement(stainValueTag());
    }
    return true;
}//end SetStainTwoRGB(std::array)

const std::array<double, 3> StainProfile::GetStainTwoRGB() const {
    //Create a zero array to return
    std::array<double, 3> out = { 0.0,0.0,0.0 };
    double r, g, b; //temp values
    //Get the first stain value element, or return false if not found
    if (m_stainTwoElement == nullptr) { return out; }
    tinyxml2::XMLElement* sVals = m_stainTwoElement->FirstChildElement(stainValueTag());
    if (sVals == nullptr) { return out; }
    while (sVals != nullptr) {
        if (sVals->Attribute(valueTypeAttribute(), "r")) {
            tinyxml2::XMLError rResult = sVals->QueryDoubleText(&r);
            if (rResult != tinyxml2::XML_SUCCESS) { return out; }
        }
        else if (sVals->Attribute(valueTypeAttribute(), "g")) {
            tinyxml2::XMLError gResult = sVals->QueryDoubleText(&g);
            if (gResult != tinyxml2::XML_SUCCESS) { return out; }
        }
        else if (sVals->Attribute(valueTypeAttribute(), "b")) {
            tinyxml2::XMLError bResult = sVals->QueryDoubleText(&b);
            if (bResult != tinyxml2::XML_SUCCESS) { return out; }
        }
        sVals = sVals->NextSiblingElement(stainValueTag());
    }
    out[0] = r;
    out[1] = g;
    out[2] = b;
    return out;
}//end GetStainTwoRGB

bool StainProfile::SetStainThreeRGB(const double &r, const double &g, const double &b) {
    //Place data into a std::array, pass to other overload of this method
    std::array<double, 3> rgb = { r,g,b };
    return SetStainThreeRGB(rgb);
}//end SetStainThreeRGB

bool StainProfile::SetStainThreeRGB(double rgb[])
{
    //Check length of rgb first
    //This is a C-style array, so get the size the old way
    int len = sizeof(rgb) / sizeof(*rgb);
    if (len != 3) {
        return false;
    }
    else {
        //Place data into a std::array, pass to other overload of this method
        std::array<double, 3> arrayRGB;
        for (int i = 0; i < 3; ++i) {
            try {
                arrayRGB.at(i) = rgb[i];
            }
            catch (const std::out_of_range& rangeerr) {
                rangeerr.what();
                //The index is out of range. Return false
                return false;
            }
        }
        return SetStainThreeRGB(arrayRGB);
    }
    return false;
}//end SetStainThreeRGB(C array)

bool StainProfile::SetStainThreeRGB(const std::array<double, 3> &rgb_in) {
    //Normalize the array before assigning
    std::array<double, 3> rgb = StainVectorMath::NormalizeArray<double, 3>(rgb_in);
    //Get the first stain value element, or return false if not found
    if (m_stainThreeElement == nullptr) { return false; }
    tinyxml2::XMLElement* sVals = m_stainThreeElement->FirstChildElement(stainValueTag());
    if (sVals == nullptr) { return false; }
    while (sVals != nullptr) {
        if (sVals->Attribute(valueTypeAttribute(), "r")) {
            sVals->SetText(rgb[0]);
        }
        else if (sVals->Attribute(valueTypeAttribute(), "g")) {
            sVals->SetText(rgb[1]);
        }
        else if (sVals->Attribute(valueTypeAttribute(), "b")) {
            sVals->SetText(rgb[2]);
        }
        sVals = sVals->NextSiblingElement(stainValueTag());
    }
    return true;
}//end SetStainThreeRGB(std::array)

const std::array<double, 3> StainProfile::GetStainThreeRGB() const {
    //Create a zero array to return
    std::array<double, 3> out = { 0.0,0.0,0.0 };
    double r, g, b; //temp values
    //Get the first stain value element, or return false if not found
    if (m_stainThreeElement == nullptr) { return out; }
    tinyxml2::XMLElement* sVals = m_stainThreeElement->FirstChildElement(stainValueTag());
    if (sVals == nullptr) { return out; }
    while (sVals != nullptr) {
        if (sVals->Attribute(valueTypeAttribute(), "r")) {
            tinyxml2::XMLError rResult = sVals->QueryDoubleText(&r);
            if (rResult != tinyxml2::XML_SUCCESS) { return out; }
        }
        else if (sVals->Attribute(valueTypeAttribute(), "g")) {
            tinyxml2::XMLError gResult = sVals->QueryDoubleText(&g);
            if (gResult != tinyxml2::XML_SUCCESS) { return out; }
        }
        else if (sVals->Attribute(valueTypeAttribute(), "b")) {
            tinyxml2::XMLError bResult = sVals->QueryDoubleText(&b);
            if (bResult != tinyxml2::XML_SUCCESS) { return out; }
        }
        sVals = sVals->NextSiblingElement(stainValueTag());
    }
    out[0] = r;
    out[1] = g;
    out[2] = b;
    return out;
}//end GetStainThreeRGB

bool StainProfile::GetProfilesAsDoubleArray(double (&profileArray)[9]) const {
    return GetProfilesAsDoubleArray(profileArray, false);
}//end GetProfileAsDoubleArray

bool StainProfile::GetNormalizedProfilesAsDoubleArray(double (&profileArray)[9]) const {
    return GetProfilesAsDoubleArray(profileArray, true);
}//end GetNormalizedProfilesAsDoubleArray

bool StainProfile::GetProfilesAsDoubleArray(double (&profileArray)[9], bool normalize) const {
    //This method fills the values of profileArray from local stain profile
    //Assigns 0.0 for elements corresponding to components beyond the number set in the profile
    int components = this->GetNumberOfStainComponents();
    if (components <= 0) {
        for (int i = 0; i < 9; i++) { profileArray[i] = 0.0; }
        return false;
    }
    std::array<double, 3> raw_rgb1 = this->GetStainOneRGB();
    std::array<double, 3> raw_rgb2 = this->GetStainTwoRGB();
    std::array<double, 3> raw_rgb3 = this->GetStainThreeRGB();

    std::array<double, 3> rgb1 = normalize ? StainVectorMath::NormalizeArray(raw_rgb1) : raw_rgb1;
    std::array<double, 3> rgb2 = normalize ? StainVectorMath::NormalizeArray(raw_rgb2) : raw_rgb2;
    std::array<double, 3> rgb3 = normalize ? StainVectorMath::NormalizeArray(raw_rgb3) : raw_rgb3;

    //Use the ternary operator for assignment of the profileArray elements
    profileArray[0] = (components > 0)  ? rgb1[0] : 0.0;
    profileArray[1] = (components > 0)  ? rgb1[1] : 0.0;
    profileArray[2] = (components > 0)  ? rgb1[2] : 0.0;

    profileArray[3] = (components > 1)  ? rgb2[0] : 0.0;
    profileArray[4] = (components > 1)  ? rgb2[1] : 0.0;
    profileArray[5] = (components > 1)  ? rgb2[2] : 0.0;

    profileArray[6] = (components > 2) ? rgb3[0] : 0.0;
    profileArray[7] = (components > 2) ? rgb3[1] : 0.0;
    profileArray[8] = (components > 2) ? rgb3[2] : 0.0;
    return true;
}//end GetProfilesAsDoubleArray

bool StainProfile::SetProfilesFromDoubleArray(double (&profileArray)[9]) {
    //Copy values to ensure they are in local scope
    double thisArray[9] = { 0.0 };
    for (int i = 0; i < 9; i++) {
        thisArray[i] = profileArray[i];
    }
    //Assume that the contents of profileArray are intentional, and directly set all values
    bool check1 = this->SetStainOneRGB(thisArray[0], thisArray[1], thisArray[2]);
    bool check2 = this->SetStainTwoRGB(thisArray[3], thisArray[4], thisArray[5]);
    bool check3 = this->SetStainThreeRGB(thisArray[6], thisArray[7], thisArray[8]);
    return (check1 && check2 && check3);
}//end SetProfileFromDoubleArray

bool StainProfile::checkFile(const std::string &fileString, const std::string &op) {
    if (fileString.empty()) {
        return false;
    }
    //Check if the file exists
    std::fstream theFile(fileString.c_str());
    bool fileExists = (bool)theFile;

    //If read, check file is readable
    //If write, check that the file can be created or written to
    if (!op.compare("r")) {
        std::ifstream inFile(fileString.c_str());
        return (inFile.good() && fileExists);
    }
    else if (!op.compare("w")) {
        std::ofstream outFile(fileString.c_str());
        bool outFileGood = outFile.good();
        bool outFileWritable = outFile.is_open();
        //If the file can be written to, return true
        if (outFileGood && outFileWritable) {
            return true;
        }
        else {
            return false;
        }
    }
    //else
    return false;
}//end checkFile

///Public write method, calls private write method
bool StainProfile::writeStainProfile(const std::string &fileString) {
    bool checkResult = this->checkFile(fileString, "w");
    if (!checkResult) {
        return false;
    }
    else {
        tinyxml2::XMLError eResult = this->writeStainProfileToXML(fileString);
        if (eResult == tinyxml2::XML_SUCCESS) {
            return true;
        }
        else {
            return false;
        }
    }
    return false;
}//end writeStainProfile

///Public read method, calls private write method
bool StainProfile::readStainProfile(const std::string &fileString) {
    bool checkResult = this->checkFile(fileString, "r");
    if (!checkResult) {
        return false;
    }
    else {
        tinyxml2::XMLError eResult = this->readStainProfileFromXMLFile(fileString);
        if (eResult == tinyxml2::XML_SUCCESS) {
            return true;
        }
        else {
            return false;
        }
    }
    return false;
}//end readStainProfile

/// Public read method, calls private write method
bool StainProfile::readStainProfile(const char *str, size_t size) {
    tinyxml2::XMLError eResult = this->readStainProfileFromXMLString(str, size);
    if (eResult == tinyxml2::XML_SUCCESS) {
      return true;
    } else {
      return false;
    }
}  // end readStainProfile

///Private write method that deals with TinyXML2
tinyxml2::XMLError StainProfile::writeStainProfileToXML(const std::string &fileString) {
    //Write it!
    tinyxml2::XMLError eResult = this->GetXMLDoc()->SaveFile(fileString.c_str());
    XMLCheckResult(eResult);
    //else
    return tinyxml2::XML_SUCCESS;
}//end writeStainProfileToXML

tinyxml2::XMLError StainProfile::readStainProfileFromXMLFile(const std::string &fileString) {
    //Clear the current structure
    bool clearSuccessful = ClearXMLDocument();
    if (!clearSuccessful) {
        return tinyxml2::XML_ERROR_PARSING_ELEMENT;
    }
    //else
    //Read it!
    tinyxml2::XMLError eResult = this->GetXMLDoc()->LoadFile(fileString.c_str());
    XMLCheckResult(eResult);

    return parseXMLDoc();
}//end readStainProfileFromXMLFile

tinyxml2::XMLError StainProfile::readStainProfileFromXMLString(const char *str, size_t size) {
  // Clear the current structure
  bool clearSuccessful = ClearXMLDocument();
  if (!clearSuccessful) {
    return tinyxml2::XML_ERROR_PARSING_ELEMENT;
  }
  // else
  // Read it!
  tinyxml2::XMLError eResult = this->GetXMLDoc()->Parse(str, size);
  XMLCheckResult(eResult);

  return parseXMLDoc();
}  // end readStainProfileFromXMLFile

tinyxml2::XMLError StainProfile::parseXMLDoc() {
    // Assign the member pointers to children in the loaded data structure
    m_rootElement = this->GetXMLDoc()->FirstChildElement(rootTag());
    if (m_rootElement == nullptr) {
        return tinyxml2::XML_ERROR_PARSING_ELEMENT;
    }
    m_componentsElement = m_rootElement->FirstChildElement(componentsTag());
    if (m_componentsElement == nullptr) {
        return tinyxml2::XML_ERROR_PARSING_ELEMENT;
    }

    tinyxml2::XMLElement *tempStain =
        m_componentsElement->FirstChildElement(stainTag());
    while (tempStain != nullptr) {
        // They must all have an index, or it should be considered an error
        int stainIndex(-1);
        tinyxml2::XMLError eResult =
            tempStain->QueryIntAttribute(indexOfStainAttribute(), &stainIndex);
        XMLCheckResult(eResult);
        // Choose which stain element to assign tempStain to
        if (stainIndex == 1) {
            m_stainOneElement = tempStain;
        } else if (stainIndex == 2) {
            m_stainTwoElement = tempStain;
        } else if (stainIndex == 3) {
            m_stainThreeElement = tempStain;
        }
        // tempStain is just a pointer to an existing element
        // Dropping it does not create a memory leak
        tempStain = tempStain->NextSiblingElement(stainTag());
    }

    m_analysisModelElement = m_rootElement->FirstChildElement(analysisModelTag());
    if (m_analysisModelElement == nullptr) {
        return tinyxml2::XML_ERROR_PARSING_ELEMENT;
    }
    //Parsing any parameters of this element is handled by Get methods

    m_algorithmElement = m_rootElement->FirstChildElement(algorithmTag());
    if (m_algorithmElement == nullptr) {
        return tinyxml2::XML_ERROR_PARSING_ELEMENT;
    }
    //Parsing any parameters of this element is handled by Get methods

    return tinyxml2::XML_SUCCESS;
}//end parseXMLDoc

const std::vector<std::string> StainProfile::GetStainAnalysisModelOptions() const {
    return m_stainAnalysisModelOptions;
}//end GetStainAnalysisModelOptions

const std::vector<std::string> StainProfile::GetStainSeparationAlgorithmOptions() const {
    return m_stainSeparationAlgorithmOptions;
}//end GetStainSeparationAlgorithmOptions

const int StainProfile::GetVectorIndexFromName(const std::string &name, const std::vector<std::string> &vec) const {
    int outIndex(-1);
    auto it = std::find(vec.begin(), vec.end(), name);
    if (it != vec.end()) {
        outIndex = static_cast<int>(std::distance(vec.begin(), it));
    }
    else {
        outIndex = -1;
    }
    return outIndex;
}//end GetVectorIndexFromName

const std::string StainProfile::GetValueFromStringVector(const int &index,
    const std::vector<std::string> &vec) const {
    //Check that the given index value is valid
    //Use the vector::at operator to do bounds checking
    std::string name;
    try {
        name = vec.at(index);
    }
    catch (const std::out_of_range& rangeerr) {
        rangeerr.what();
        //The index is out of range. Return empty string
        return std::string();
    }
    //name found. Return it.
    return name;
}//end GetValueFromStringVector

const std::string StainProfile::GetStainAnalysisModelName(const int &index) const {
    return GetValueFromStringVector(index, m_stainAnalysisModelOptions);
}//end GetStainAnalysisModelName

const std::string StainProfile::GetStainSeparationAlgorithmName(const int &index) const {
    return GetValueFromStringVector(index, m_stainSeparationAlgorithmOptions);
}//end GetStainSeparationAlgorithmName

bool StainProfile::BuildXMLDocument() {
    //Instantiate the XMLDocument as shared_ptr
    m_xmlDoc = std::make_shared<tinyxml2::XMLDocument>();
    //Create a root element and assign it as a child of the XML document
    m_rootElement = m_xmlDoc->NewElement(rootTag());
    m_rootElement->SetAttribute(nameOfStainProfileAttribute(), "");
    m_xmlDoc->InsertFirstChild(m_rootElement);

    //Build the next levels of the document structure
    m_componentsElement = m_xmlDoc->NewElement(componentsTag());
    m_componentsElement->SetAttribute(numberOfStainsAttribute(), "0"); //default value 0
    m_rootElement->InsertEndChild(m_componentsElement);

    //Build stain one
    m_stainOneElement = m_xmlDoc->NewElement(stainTag());
    m_stainOneElement->SetAttribute(indexOfStainAttribute(), 1);
    m_stainOneElement->SetAttribute(nameOfStainAttribute(), "");
    m_componentsElement->InsertEndChild(m_stainOneElement);
    //Add the three stain values to stain element one
    //r
    tinyxml2::XMLElement* rValOne = m_xmlDoc->NewElement(stainValueTag());
    rValOne->SetAttribute(valueTypeAttribute(), "r");
    rValOne->SetText(0);
    m_stainOneElement->InsertEndChild(rValOne);
    //g
    tinyxml2::XMLElement* gValOne = m_xmlDoc->NewElement(stainValueTag());
    gValOne->SetAttribute(valueTypeAttribute(), "g");
    gValOne->SetText(0);
    m_stainOneElement->InsertEndChild(gValOne);
    //b
    tinyxml2::XMLElement* bValOne = m_xmlDoc->NewElement(stainValueTag());
    bValOne->SetAttribute(valueTypeAttribute(), "b");
    bValOne->SetText(0);
    m_stainOneElement->InsertEndChild(bValOne);

    //Build stain two
    m_stainTwoElement = m_xmlDoc->NewElement(stainTag());
    m_stainTwoElement->SetAttribute(indexOfStainAttribute(), 2);
    m_stainTwoElement->SetAttribute(nameOfStainAttribute(), "");
    m_componentsElement->InsertEndChild(m_stainTwoElement);
    //Add the three stain values to stain element two
    //r
    tinyxml2::XMLElement* rValTwo = m_xmlDoc->NewElement(stainValueTag());
    rValTwo->SetAttribute(valueTypeAttribute(), "r");
    rValTwo->SetText(0);
    m_stainTwoElement->InsertEndChild(rValTwo);
    //g
    tinyxml2::XMLElement* gValTwo = m_xmlDoc->NewElement(stainValueTag());
    gValTwo->SetAttribute(valueTypeAttribute(), "g");
    gValTwo->SetText(0);
    m_stainTwoElement->InsertEndChild(gValTwo);
    //b
    tinyxml2::XMLElement* bValTwo = m_xmlDoc->NewElement(stainValueTag());
    bValTwo->SetAttribute(valueTypeAttribute(), "b");
    bValTwo->SetText(0);
    m_stainTwoElement->InsertEndChild(bValTwo);

    //Build stain three
    m_stainThreeElement = m_xmlDoc->NewElement(stainTag());
    m_stainThreeElement->SetAttribute(indexOfStainAttribute(), 3);
    m_stainThreeElement->SetAttribute(nameOfStainAttribute(), "");
    m_componentsElement->InsertEndChild(m_stainThreeElement);
    //Add the three stain values to stain element three
    //r
    tinyxml2::XMLElement* rValThree = m_xmlDoc->NewElement(stainValueTag());
    rValThree->SetAttribute(valueTypeAttribute(), "r");
    rValThree->SetText(0);
    m_stainThreeElement->InsertEndChild(rValThree);
    //g
    tinyxml2::XMLElement* gValThree = m_xmlDoc->NewElement(stainValueTag());
    gValThree->SetAttribute(valueTypeAttribute(), "g");
    gValThree->SetText(0);
    m_stainThreeElement->InsertEndChild(gValThree);
    //b
    tinyxml2::XMLElement* bValThree = m_xmlDoc->NewElement(stainValueTag());
    bValThree->SetAttribute(valueTypeAttribute(), "b");
    bValThree->SetText(0);
    m_stainThreeElement->InsertEndChild(bValThree);

    //analysis model tag
    m_analysisModelElement = m_xmlDoc->NewElement(analysisModelTag());
    m_rootElement->InsertEndChild(m_analysisModelElement);

    //algorithm tag, which can contain parameter values needed by the algorithm
    m_algorithmElement = m_xmlDoc->NewElement(algorithmTag());
    m_rootElement->InsertEndChild(m_algorithmElement);

    return true;
}//end BuildXMLDocument

bool StainProfile::CheckXMLDocument() {
    //Check if the basic structure of this XMLDocument is complete
    //Get the pointer to the xmldocument
    std::shared_ptr<tinyxml2::XMLDocument> docPtr = this->GetXMLDoc();
    if (docPtr == nullptr) {
        return false;
    }
    if (docPtr->NoChildren()) {
        return false;
    }
    //else

    //Check that all of the 'Get' operations from this class will be successful
    //TODO: additional checks
    try {
        auto a = this->GetNameOfStainProfile();
        auto b = this->GetNumberOfStainComponents();
        auto c = this->GetNameOfStainOne();
        auto d = this->GetNameOfStainTwo();
        auto e = this->GetNameOfStainThree();
        auto f = this->GetNameOfStainAnalysisModel();
        auto g = this->GetNameOfStainSeparationAlgorithm();
        auto h = this->GetStainOneRGB();
        auto i = this->GetStainTwoRGB();
        auto j = this->GetStainThreeRGB();
        double k[9] = { 0.0 };
        bool l = this->GetProfilesAsDoubleArray(k, true);
        bool m = this->GetProfilesAsDoubleArray(k, false);
        if (!l || !m) {
            return false;
        }
    }
    catch (...) {
        return false;
    }
    //else
    return true;
}//end CheckXMLDocument

bool StainProfile::ClearXMLDocument() {
    //Set all of the values in the member XML document to their null values
    std::shared_ptr<tinyxml2::XMLDocument> docPtr = this->GetXMLDoc();
    //else
    //Use the member methods
    this->SetNameOfStainProfile("");
    this->SetNumberOfStainComponents(0);
    this->SetNameOfStainOne("");
    this->SetNameOfStainTwo("");
    this->SetNameOfStainThree("");
    this->SetNameOfStainAnalysisModel("");
    this->SetNameOfStainSeparationAlgorithm(""); 

    //Clear the stain vector values using the specific method in this class
    this->ClearStainVectorValues();
    //Clear the analysis model parameters (if any)
    this->ClearAllAnalysisModelParameters();
    //Clear the separation algorithm parameters (if any)
    this->ClearAllSeparationAlgorithmParameters();

    //Success
    return true;
}//end ClearXMLDocument

bool StainProfile::ClearStainVectorValues() {
    //Set the stain vector values in the member XML document to 0
    std::shared_ptr<tinyxml2::XMLDocument> docPtr = this->GetXMLDoc();
    this->SetStainOneRGB(0, 0, 0);
    this->SetStainTwoRGB(0, 0, 0);
    this->SetStainThreeRGB(0, 0, 0);

    //Success
    return true;
}//end ClearStainVectorValues

bool StainProfile::ClearChildren(tinyxml2::XMLElement* el, const std::string &tag /*= std::string()*/) {
    bool success = true;
    bool errorVal = false;
    if (el != nullptr) {
        if (!el->NoChildren()) {
            if (tag.empty()) {
                //If no tag was specified, delete all children
                el->DeleteChildren();
                return success;
            }
            else {
                while (RemoveSingleChild(el, tag)) {
                    ;//as long as RemoveSingleChild returns true, repeat running it
                }
                return success;
            }
        }//else
        return errorVal;
    }//else
    return errorVal;
}//end ClearChildren

bool StainProfile::ClearAllAnalysisModelParameters() {
    std::string tag = std::string(parameterTag());
    return ClearChildren(m_analysisModelElement, tag);
}//end ClearAnalysisModelParameters

bool StainProfile::ClearAllSeparationAlgorithmParameters() {
    std::string tag = std::string(parameterTag());
    return ClearChildren(m_algorithmElement, tag);
}//end ClearSeparationAlgorithmParameters

bool StainProfile::RemoveSingleChild(tinyxml2::XMLElement* el, const std::string &tag /*= std::string()*/,
    const std::string &att /*= std::string()*/, const std::string &val /*= std::string()*/) {
    bool success = true;
    bool errorVal = false;
    tinyxml2::XMLElement* child;
    if (el != nullptr) {
        if (!el->NoChildren()) {
            //Get a reference to the child to be deleted, if it can be found
            if (tag.empty()) {
                //If tag is an empty string, the child to be deleted is the first one,
                //without regard to its xml tag
                child = el->FirstChildElement();
                //If there is no valid child, return errorVal
                if (child == nullptr) { return errorVal; }
            }
            else {
                //If no attribute is specified, choose the first child with the given tag
                if (att.empty()) {
                    child = el->FirstChildElement(tag.c_str());
                }
                else {
                    tinyxml2::XMLElement* tempChild = el->FirstChildElement(tag.c_str());
                    //loop through all child elements with this tag
                    bool attFound = false;
                    while (tempChild != nullptr) {
                        if (val.empty()) {
                            //If the parameter has the given attribute set to anything,
                            //set attFound to true
                            const char* tempAtt = tempChild->Attribute(att.c_str());
                            if (tempAtt != nullptr) { attFound = true; }
                        }
                        else {
                            //If there is a requested value for the attribute, test with it
                            const char* tempAtt = tempChild->Attribute(att.c_str(), val.c_str());
                            if (tempAtt != nullptr) { attFound = true; }
                        }
                        if (attFound) {
                            child = tempChild;
                            break;
                        }
                        //else
                        tempChild = tempChild->NextSiblingElement(tag.c_str());
                    }//end while
                }
            }

            //If a child to be deleted is found, delete it, otherwise return errorVal
            if (child != nullptr) {
                el->DeleteChild(child);
                return success;
            }
            //else
            return errorVal;
        }//else
        return errorVal;
    }//else
    return errorVal;
}//end RemoveSingleChild

bool StainProfile::RemoveAnalysisModelParameter(const std::string &pType) {
    std::string paramTag = std::string(parameterTag());
    std::string paramAtt = std::string(parameterTypeAttribute());
    return RemoveSingleChild(m_analysisModelElement, paramTag, paramAtt, pType);
}//end RemoveAnalysisModelParameter

bool StainProfile::RemoveSeparationAlgorithmParameter(const std::string &pType) {
    std::string paramTag = std::string(parameterTag());
    std::string paramAtt = std::string(parameterTypeAttribute());
    return RemoveSingleChild(m_algorithmElement, paramTag, paramAtt, pType);
}//end RemoveSeparationAlgorithmParameter

const std::map<std::string, std::string> StainProfile::GetAllParameters(tinyxml2::XMLElement* el) const {
    auto errorVal = std::map<std::string, std::string>();
    auto theMap = std::map<std::string, std::string>();
    if (el == nullptr) { return errorVal; } //if nothing can be set, return an empty map
    if (m_xmlDoc == nullptr) { return errorVal; } //if there is no xml document, return errorVal
    //This gets the value of a parameter with param-type given by the type argument, if defined
    const char* paramTag = parameterTag();
    const char* paramAtt = parameterTypeAttribute();

    //Get the first child element, or return empty string if not found
    tinyxml2::XMLElement* tempParam = el->FirstChildElement(paramTag);
    //Return errorVal if there are no parameter tags
    if (tempParam == nullptr) { return errorVal; }
    while (tempParam != nullptr) {
        const char* tempAtt = tempParam->Attribute(paramAtt);
        const char* tempVal = tempParam->GetText();
        //If both values exist, create a pair and append to the map
        if ((tempAtt != nullptr) && (tempVal != nullptr)) {
            std::pair<std::string, std::string> tempPair
                = std::pair<std::string, std::string>(std::string(tempAtt), std::string(tempVal));
            theMap.insert(tempPair);
        }
        //repeat loop until there are no more parameters
        tempParam = tempParam->NextSiblingElement(paramTag);
    }
    return theMap;
}//end GetAllParameters

bool StainProfile::SetAllParameters(tinyxml2::XMLElement* el, const std::map<std::string, std::string>& p) {
    if (p.empty()) { return false; }     //if nothing can be set, return false
    if (el == nullptr) { return false; } //if nothing can be set, return false
    //Clear all children of the element node before assignment
    this->ClearChildren(el);

    //Iterate over map, create and assign all parameter element nodes
    for (auto it = p.begin(); it != p.end(); ++it) {
        std::string key = it->first;
        std::string val = it->second;
        bool checkVal = this->SetSingleParameter(el, key, val);
        if (!checkVal) { return false; }
    }
    return true;
}//end SetAllParameters

const std::map<std::string, std::string> StainProfile::GetAllAnalysisModelParameters() const {
    return GetAllParameters(m_analysisModelElement);
}//end GetAnalysisModelParameters

bool StainProfile::SetAllAnalysisModelParameters(const std::map<std::string, std::string>& p) {
    return SetAllParameters(m_analysisModelElement, p);
}//end SetAnalysisModelParameters

const std::map<std::string, std::string> StainProfile::GetAllSeparationAlgorithmParameters() const {
    return GetAllParameters(m_algorithmElement);
}//end GetSeparationAlgorithmParameters

bool StainProfile::SetAllSeparationAlgorithmParameters(const std::map<std::string, std::string>& p) {
    return SetAllParameters(m_algorithmElement, p);
}//end SetSeparationAlgorithmParameters

const std::string StainProfile::GetSingleParameter(tinyxml2::XMLElement* el, const std::string &type) const {
    std::string errorVal = std::string(); //an empty string
    std::string outString = std::string();
    if (el == nullptr) { return errorVal; } //if there is no element to get from, return errorVal
    if (m_xmlDoc == nullptr) { return errorVal; } //if there is no xml document, return errorVal
    if (type.empty()) { return errorVal; } //if no type is specified, return errorVal

    //This gets the value of a parameter with param-type given by the type argument, if defined
    const char* paramTag = parameterTag();
    const char* paramAtt = parameterTypeAttribute();

    //Get the first child element, or return empty string if not found
    tinyxml2::XMLElement* tempParam = el->FirstChildElement(paramTag);
    //Return errorVal if there are no parameter tags
    if (tempParam == nullptr) { return errorVal; }
    while (tempParam != nullptr) {
        if (tempParam->Attribute(paramAtt, type.c_str())) {
            const char *buffer = tempParam->Value();
            if (buffer == nullptr) {
                //The parameter of the right type was found, but has no value set
                outString = errorVal;
            }
            else {
                outString = std::string(buffer);
            }
            break;
        }
        //else repeat loop
        tempParam = tempParam->NextSiblingElement(paramTag);
    }
    return outString;
}//end GetSingleParameter

bool StainProfile::SetSingleParameter(tinyxml2::XMLElement* el, const std::string &type, const std::string &val) {
    if (type.empty()) { return false; }   //if there is no type, nothing can be set; return false
    if (el == nullptr) { return false; } //if nothing can be set, return false
    if (m_xmlDoc == nullptr) { return false; } //if there is no xml document, return false
    //This sets the value of a parameter with param-type given by the type argument
    const char* paramTag = parameterTag();
    const char* paramAtt = parameterTypeAttribute();
    bool attFound = false;
    //Get the first child element, if none while loop will not be entered
    tinyxml2::XMLElement* tempParam = el->FirstChildElement(paramTag);
    //If at least one parameter exists, loop through to find parameter of given type
    while (tempParam != nullptr) {
        const char* tempAtt = tempParam->Attribute(paramAtt, type.c_str());
        if (tempAtt != nullptr) {
            tempParam->SetText(val.c_str());
            attFound = true;
        }
        //else
        //move to next parameter, prepare to repeat
        tempParam = tempParam->NextSiblingElement(paramTag);
    }
    //if while loop was entered and attribute of given type was not found, create and add it
    if(!attFound) {
        //Create a new parameter, add to XML document model
        tinyxml2::XMLElement* newParam = m_xmlDoc->NewElement(parameterTag());
        newParam->SetAttribute(parameterTypeAttribute(), type.c_str());
        newParam->SetText(val.c_str());
        el->InsertEndChild(newParam);
    }
    return true;
}//end SetSingleParameter

const std::string StainProfile::GetSingleAnalysisModelParameter(const std::string &type) const {
    return GetSingleParameter(m_analysisModelElement, type);
}//end GetSingleAnalysisModelParameter

bool StainProfile::SetSingleAnalysisModelParameter(const std::string &type, const std::string &val) {
    return SetSingleParameter(m_analysisModelElement, type, val);
}//end SetSingleAnalysisModelParameter

const std::string StainProfile::GetSingleSeparationAlgorithmParameter(const std::string &type) const {
    return GetSingleParameter(m_algorithmElement, type);
}//end GetSingleSeparationAlgorithmParameter

bool StainProfile::SetSingleSeparationAlgorithmParameter(const std::string &type, const std::string &val) {
    return SetSingleParameter(m_algorithmElement, type, val);
}//end SetSingleSeparationAlgorithmParameter

const long int StainProfile::GetSeparationAlgorithmNumPixelsParameter() const {
    std::string type = this->pTypeNumPixels();
    std::string oss = this->GetSingleSeparationAlgorithmParameter(type);
    long int outVal(-1); //Set error value here
    //Convert oss to long int, if possible. Catch all exceptions (invalid_argument and out_of_range)
    try {
        outVal = std::stol(oss);
    }
    catch (...) {
        //No additional actions required
    }
    return outVal;
}//end GetSeparationAlgorithmNumPixelsParameter

bool StainProfile::SetSeparationAlgorithmNumPixelsParameter(const long int& p) {
    std::string type = this->pTypeNumPixels();
    std::stringstream ss;
    ss << p; // << std::endl;
    return this->SetSingleSeparationAlgorithmParameter(type, ss.str());
}//end SetSeparationAlgorithmNumPixelsParameter

const double StainProfile::GetSeparationAlgorithmThresholdParameter() const {
    std::string type = this->pTypeThreshold();
    std::string oss = this->GetSingleSeparationAlgorithmParameter(type);
    double outVal(-1.); //Set error value here
    //Convert oss to double, if possible. Catch all exceptions (invalid_argument and out_of_range)
    try {
        outVal = std::stod(oss);
    }
    catch (...) {
        //No additional actions required
    }
    return outVal;
}//end GetSeparationAlgorithmThresholdParameter

bool StainProfile::SetSeparationAlgorithmThresholdParameter(const double& p) {
    std::string type = this->pTypeThreshold();
    std::stringstream ss;
    ss << p; // << std::endl;
    return this->SetSingleSeparationAlgorithmParameter(type, ss.str());
}//end SetSeparationAlgorithmThresholdParameter

const double StainProfile::GetSeparationAlgorithmPercentileParameter() const {
    std::string type = this->pTypePercentile();
    std::string oss = this->GetSingleSeparationAlgorithmParameter(type);
    double outVal(-1.); //Set error value here
    //Convert oss to double, if possible. Catch all exceptions (invalid_argument and out_of_range)
    try {
        outVal = std::stod(oss);
    }
    catch (...) {
        //No additional actions required
    }
    return outVal;
}//end GetSeparationAlgorithmPercentileParameter

bool StainProfile::SetSeparationAlgorithmPercentileParameter(const double& p) {
    std::string type = this->pTypePercentile();
    std::stringstream ss;
    ss << p; // << std::endl;
    return this->SetSingleSeparationAlgorithmParameter(type, ss.str());
}//end SetSeparationAlgorithmPercentileParameter

const int StainProfile::GetSeparationAlgorithmHistogramBinsParameter() const {
    std::string type = this->pTypeHistoBins();
    std::string oss = this->GetSingleSeparationAlgorithmParameter(type);
    int outVal(-1); //Set error value here
    //Convert oss to int, if possible. Catch all exceptions (invalid_argument and out_of_range)
    try {
        outVal = std::stoi(oss);
    }
    catch (...) {
        //No additional actions required
    }
    return outVal;
}//end GetSeparationAlgorithmHistogramBinsParameter

bool StainProfile::SetSeparationAlgorithmHistogramBinsParameter(const int& p) {
    std::string type = this->pTypeHistoBins();
    std::stringstream ss;
    ss << p; // << std::endl;
    return this->SetSingleSeparationAlgorithmParameter(type, ss.str());
}//end SetSeparationAlgorithmHistogramBinsParameter

