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

#ifndef SEDEEN_SRC_PLUGINS_STAINANALYSIS_STAINPROFILE_H
#define SEDEEN_SRC_PLUGINS_STAINANALYSIS_STAINPROFILE_H

#include <string>
#include <vector>
#include <array>
#include <memory>

#include "StainVectorMath.h"

//TinyXML2 system for XML handling
#include <tinyxml2.h>
//Macro to exit from a method/function if an XML operation was unsuccessful
#ifndef XMLCheckResult
    #define XMLCheckResult(a_eResult) if (a_eResult != tinyxml2::XML_SUCCESS) { return a_eResult; }
#endif

///Class to read and write stain vector profile data to XML files. 
///Avoids invoking the Sedeen SDK. Uses TinyXML2 for XML file interactions.
///Stores data in the member XMLDocument structure
class StainProfile 
{
public:
    ///No-parameter constructor (builds XMLDocument structure)
    StainProfile();
    ///Copy constructor
    StainProfile(StainProfile &s);
    ///virtual destructor
    virtual ~StainProfile();

    ///Get/Set the name of the stain profile
    bool SetNameOfStainProfile(const std::string &);
    ///Get/Set the name of the stain profile
    const std::string GetNameOfStainProfile();

    ///Get/Set the number of stain components in the profile
    bool SetNumberOfStainComponents(const int &);
    ///Get/Set the number of stain components in the profile
    const int GetNumberOfStainComponents();

    ///Get/Set the name of stain one
    bool SetNameOfStainOne(const std::string &);
    ///Get/Set the name of stain one
    const std::string GetNameOfStainOne();

    ///Get/Set the name of stain two
    bool SetNameOfStainTwo(const std::string &);
    ///Get/Set the name of stain two
    const std::string GetNameOfStainTwo();

    ///Get/Set the name of stain three
    bool SetNameOfStainThree(const std::string &);
    ///Get/Set the name of stain three
    const std::string GetNameOfStainThree();

    ///Get/Set the name of the stain separation algorithm currently selected
    bool SetNameOfStainSeparationAlgorithm(const std::string &);
    ///Get/Set the name of the stain separation algorithm currently selected
    const std::string GetNameOfStainSeparationAlgorithm();

    ///Set the RGB values for stain one as three doubles
    bool SetStainOneRGB(double, double, double);
    ///Overload: Set the RGB values for stain one as a three-element C-style array
    bool SetStainOneRGB(double[]);
    ///Overload: Set the RGB values for stain one as a C++11 std::array of three doubles
    bool SetStainOneRGB(const std::array<double, 3> &);
    ///Get the RGB values for stain one as a C++11 std::array of three doubles
    const std::array<double, 3> GetStainOneRGB();

    ///Set the RGB values for stain two as three doubles
    bool SetStainTwoRGB(double, double, double);
    ///Overload: Set the RGB values for stain two as a three-element C-style array
    bool SetStainTwoRGB(double[]);
    ///Overload: Set the RGB values for stain two as a C++11 std::array of three doubles
    bool SetStainTwoRGB(const std::array<double, 3> &);
    ///Get the RGB values for stain two as a C++11 std::array of three doubles
    const std::array<double, 3> GetStainTwoRGB();

    ///Set the RGB values for stain three as three doubles
    bool SetStainThreeRGB(double, double, double);
    ///Overload: Set the RGB values for stain three as a three-element C-style array
    bool SetStainThreeRGB(double[]);
    ///Overload: Set the RGB values for stain three as a C++11 std::array of three doubles
    bool SetStainThreeRGB(const std::array<double, 3> &);
    ///Get the RGB values for stain three as a C++11 std::array of three doubles
    const std::array<double, 3> GetStainThreeRGB();

    ///Public write method - returns true/false for success/failure
    bool writeStainProfile(std::string);
    ///Public read method - checks file existence first, returns true/false for success/failure
    bool readStainProfile(std::string);
    /// Public read method - reads profile from embedded resource
    bool readStainProfile(const char *, size_t);

    ///Request the list of possible stain separation algorithm names from the class, static defined in constructor
    const std::vector<std::string> GetStainSeparationAlgorithmOptions();

    ///Request an element of the vector of stain separation algorithms. Returns "" on error.
    const std::string GetStainSeparationAlgorithmName(int);

    ///Get the raw (no normalization applied) stain vector profiles and assign to a 9-element double array
    bool GetProfilesAsDoubleArray(double[9]);
    ///Normalize each vector to 1, then assign the stain vector profiles to a 9-element double array
    bool GetNormalizedProfilesAsDoubleArray(double[9]);
    ///Get the stain vector profiles, normalized or not depending on second parameter, and assign to a 9-element double array
    bool GetProfilesAsDoubleArray(double[9], bool);
    ///Set the values of the stain vector profiles from a 9-element double array
    bool SetProfilesFromDoubleArray(double[9]);

    //Check if the file exists, and accessible for reading or writing, depending on the second argument
    static bool checkFile(std::string, std::string);

    ///Check if the basic structure of the stain profile has been assembled
    inline bool CheckProfile() { return CheckXMLDocument(); }
    ///Clear the entire contents of the stain profile
    inline bool ClearProfile() { return ClearXMLDocument(); }
    ///Clear the stain vector values only (no changes to text fields)
    bool ClearStainVectorValues();

public:
    //XML tag strings
    //The root tag and one attribute
    static inline const char* rootTag() { return "stain-profile"; }
    static inline const char* nameOfStainProfileAttribute() { return "profile-name"; }
    //The components tag and one attribute
    static inline const char* componentsTag() { return "components"; }
    static inline const char* numberOfStainsAttribute() { return "numstains"; }
    //The stain tag and two attributes
    static inline const char* stainTag() { return "stain"; }
    static inline const char* indexOfStainAttribute() { return "index"; }
    static inline const char* nameOfStainAttribute() { return "stain-name"; }
    //The stain value tag and one attribute
    static inline const char* stainValueTag() { return "stain-value"; }
    static inline const char* valueTypeAttribute() { return "value-type"; }
    //The algorithm tag and one attribute
    static inline const char* algorithmTag() { return "algorithm"; }
    static inline const char* algorithmNameAttribute() { return "alg-name"; }
    //The parameter tag and one attribute
    static inline const char* parameterTag() { return "parameter"; }
    static inline const char* parameterTypeAttribute() { return "param-type"; }

private:
    ///Build the XMLDocument data structure
    bool BuildXMLDocument();
    ///Check if the basic structure of the member XMLDocument has been assembled
    bool CheckXMLDocument();
    ///Clear the entire contents of the member XMLDocument
    bool ClearXMLDocument();

    ///If file is able to be opened, write the current stain profile to file as XML
    tinyxml2::XMLError writeStainProfileToXML(std::string);
    ///If file is able to be opened, read from an XML file and fill variables in this class
    tinyxml2::XMLError readStainProfileFromXMLFile(std::string);
    ///Read a stain profile from a string
    tinyxml2::XMLError readStainProfileFromXMLString(const char *, size_t);
    //Parse XML document
    tinyxml2::XMLError parseXMLDoc();

private:
    ///Store the list of possible stain separation algorithm names here
	const std::vector<std::string> m_stainSeparationAlgorithmOptions;

    ///An XML document associated with this class: note that elements can't be smartpointers
    std::shared_ptr<tinyxml2::XMLDocument> m_xmlDoc;
    ///Get pointer to the member m_xmlDoc
    inline std::shared_ptr<tinyxml2::XMLDocument> GetXMLDoc() { return m_xmlDoc; }
    ///root element: cannot be a smartpointer because destructor is private and XMLDocument must handle memory management
    tinyxml2::XMLElement* m_rootElement;
    ///components element
    tinyxml2::XMLElement* m_componentsElement;
    ///stain elements
    tinyxml2::XMLElement* m_stainOneElement;
    tinyxml2::XMLElement* m_stainTwoElement;
    tinyxml2::XMLElement* m_stainThreeElement;
    ///algorithm element
    tinyxml2::XMLElement* m_algorithmElement;
    ///parameter element
    tinyxml2::XMLElement* m_parameterElement;

};

#endif