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

#ifndef SEDEEN_SRC_PLUGINS_STAINANALYSIS_STAINPROFILE_H
#define SEDEEN_SRC_PLUGINS_STAINANALYSIS_STAINPROFILE_H

#include <string>
#include <vector>
#include <array>
#include <memory>
#include <map>

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
    ///destructor
    virtual ~StainProfile();

    ///Get/Set the name of the stain profile
    bool SetNameOfStainProfile(const std::string &);
    ///Get/Set the name of the stain profile
    const std::string GetNameOfStainProfile() const;

    ///Get/Set the number of stain components in the profile
    bool SetNumberOfStainComponents(const int &);
    ///Get/Set the number of stain components in the profile
    const int GetNumberOfStainComponents() const;

    ///Get/Set the name of stain one
    bool SetNameOfStainOne(const std::string &);
    ///Get/Set the name of stain one
    const std::string GetNameOfStainOne() const;

    ///Get/Set the name of stain two
    bool SetNameOfStainTwo(const std::string &);
    ///Get/Set the name of stain two
    const std::string GetNameOfStainTwo() const;

    ///Get/Set the name of stain three
    bool SetNameOfStainThree(const std::string &);
    ///Get/Set the name of stain three
    const std::string GetNameOfStainThree() const;

    ///Get/Set the name of the stain analysis model currently selected
    bool SetNameOfStainAnalysisModel(const std::string &);
    ///Get/Set the name of the stain analysis model currently selected
    const std::string GetNameOfStainAnalysisModel() const;

    ///Get/Set the name of the stain separation method currently selected
    bool SetNameOfStainSeparationAlgorithm(const std::string &);
    ///Get/Set the name of the stain separation method currently selected
    const std::string GetNameOfStainSeparationAlgorithm() const;

    ///Set the RGB values for stain one as three doubles
    bool SetStainOneRGB(const double&, const double&, const double&);
    ///Overload: Set the RGB values for stain one as a three-element C-style array
    bool SetStainOneRGB(double[]);
    ///Overload: Set the RGB values for stain one as a C++11 std::array of three doubles
    bool SetStainOneRGB(const std::array<double, 3> &);
    ///Get the RGB values for stain one as a C++11 std::array of three doubles
    const std::array<double, 3> GetStainOneRGB() const;

    ///Set the RGB values for stain two as three doubles
    bool SetStainTwoRGB(const double&, const double&, const double&);
    ///Overload: Set the RGB values for stain two as a three-element C-style array
    bool SetStainTwoRGB(double[]);
    ///Overload: Set the RGB values for stain two as a C++11 std::array of three doubles
    bool SetStainTwoRGB(const std::array<double, 3> &);
    ///Get the RGB values for stain two as a C++11 std::array of three doubles
    const std::array<double, 3> GetStainTwoRGB() const;

    ///Set the RGB values for stain three as three doubles
    bool SetStainThreeRGB(const double&, const double&, const double&);
    ///Overload: Set the RGB values for stain three as a three-element C-style array
    bool SetStainThreeRGB(double[]);
    ///Overload: Set the RGB values for stain three as a C++11 std::array of three doubles
    bool SetStainThreeRGB(const std::array<double, 3> &);
    ///Get the RGB values for stain three as a C++11 std::array of three doubles
    const std::array<double, 3> GetStainThreeRGB() const;

public:
    ///Public write method - returns true/false for success/failure
    bool writeStainProfile(const std::string &);
    ///Public read method - checks file existence first, returns true/false for success/failure
    bool readStainProfile(const std::string &);
    /// Public read method - reads profile from embedded resource
    bool readStainProfile(const char *, size_t);

    ///General method for retrieving the index position of a given string within a vector of strings, returns -1 if not found
    const int GetVectorIndexFromName(const std::string &name, const std::vector<std::string> &vec) const;

    ///General method for retrieving a name from a given vector at a given index, protected from range errors.
    const std::string StainProfile::GetValueFromStringVector(const int &index, const std::vector<std::string> &vec) const;

    ///Request the list of possible stain analysis model names from the class (presently one option)
    const std::vector<std::string> GetStainAnalysisModelOptions() const;

    ///Request an element of the vector of stain analysis models. Returns std::string() on error.
    const std::string GetStainAnalysisModelName(const int &index) const;

    ///Request the list of possible stain separation algorithm names from the class, static defined in constructor
    const std::vector<std::string> GetStainSeparationAlgorithmOptions() const;

    ///Request an element of the vector of stain separation algorithms. Returns std::string() on error.
    const std::string GetStainSeparationAlgorithmName(const int &index) const;

    ///Get the raw (no normalization applied) stain vector profiles and assign to a 9-element double array
    bool GetProfilesAsDoubleArray(double (&profileArray)[9]) const;
    ///Normalize each vector to 1, then assign the stain vector profiles to a 9-element double array
    bool GetNormalizedProfilesAsDoubleArray(double (&profileArray)[9]) const;
    ///Get the stain vector profiles, normalized or not depending on second parameter, and assign to a 9-element double array
    bool GetProfilesAsDoubleArray(double (&profileArray)[9], bool) const;
    ///Set the values of the stain vector profiles from a 9-element double array
    bool SetProfilesFromDoubleArray(double (&profileArray)[9]);

    //Check if the file exists and accessible for reading or writing, or that the directory to write to exists
    static bool checkFile(const std::string &, const std::string &);

    ///Check if the basic structure of the stain profile has been assembled
    inline bool CheckProfile() { return CheckXMLDocument(); }
    ///Clear the entire contents of the stain profile
    inline bool ClearProfile() { return ClearXMLDocument(); }
    ///Clear the stain vector values only (no changes to text fields)
    bool ClearStainVectorValues();

    ///Remove all of the child nodes of an element, either all or of a specified tag
    bool ClearChildren(tinyxml2::XMLElement* el, const std::string &tag = std::string());
    ///Clear the analysis model parameters (remove nodes)
    bool ClearAllAnalysisModelParameters();
    ///Clear the separation algorithm parameters (remove nodes)
    bool ClearAllSeparationAlgorithmParameters();

    ///Remove the first child node of an element with a optional tag, optional attribute, and optional attribute value
    bool RemoveSingleChild(tinyxml2::XMLElement* el, const std::string &tag = std::string(), 
        const std::string &att = std::string(), const std::string &val = std::string());

    ///Remove a single analysis model parameter node with the given param-type value
    bool RemoveAnalysisModelParameter(const std::string &pType);
    ///Remove a single separation algorithm parameter node with the given param-type value
    bool RemoveSeparationAlgorithmParameter(const std::string &pType);

public:
    ///Get all of the parameters for a given element, return empty map on error or no parameters
    const std::map<std::string, std::string> GetAllParameters(tinyxml2::XMLElement* el) const;
    ///Set all of the parameters for a given element, clear and replace all parameters
    bool SetAllParameters(tinyxml2::XMLElement* el, const std::map<std::string, std::string>& p);

    ///Get all of the stain analysis model parameters
    const std::map<std::string, std::string> GetAllAnalysisModelParameters() const;
    ///Set all of the stain analysis model parameters, clear any current values
    bool SetAllAnalysisModelParameters(const std::map<std::string, std::string>& p);

    ///Get all of the separation algorithm parameters
    const std::map<std::string, std::string> GetAllSeparationAlgorithmParameters() const;
    ///Set all of the separation algorithm parameters, clear any current values
    bool SetAllSeparationAlgorithmParameters(const std::map<std::string, std::string>& p);

    ///Get a single parameter with a given param-type from the given element, return empty string if not found
    const std::string GetSingleParameter(tinyxml2::XMLElement* el, const std::string &type) const;
    ///Set a single parameter with a given param-type as a child of the given element, create if it does not exist
    bool SetSingleParameter(tinyxml2::XMLElement* el, const std::string &type, const std::string &val);

    ///Get a single stain analysis model parameter, return empty string if not found
    const std::string GetSingleAnalysisModelParameter(const std::string &type) const;
    ///Set a single stain analysis model parameter, create if it does not exist
    bool SetSingleAnalysisModelParameter(const std::string &type, const std::string &val);

    ///Get a single stain separation algorithm parameter, return empty string if not found
    const std::string GetSingleSeparationAlgorithmParameter(const std::string &type) const;
    ///Set a single separation algorithm parameter, create if it does not exist
    bool SetSingleSeparationAlgorithmParameter(const std::string &type, const std::string &val);

    ///Set/Get the SeparationAlgorithm NumPixels parameter
    const long int GetSeparationAlgorithmNumPixelsParameter() const;
    ///Set/Get the SeparationAlgorithm NumPixels parameter
    bool SetSeparationAlgorithmNumPixelsParameter(const long int& p);

    ///Set/Get the SeparationAlgorithm Threshold parameter
    const double GetSeparationAlgorithmThresholdParameter() const;
    ///Set/Get the SeparationAlgorithm Threshold parameter
    bool SetSeparationAlgorithmThresholdParameter(const double & p);

    ///Set/Get the SeparationAlgorithm Percentile parameter
    const double GetSeparationAlgorithmPercentileParameter() const;
    ///Set/Get the SeparationAlgorithm Percentile parameter
    bool SetSeparationAlgorithmPercentileParameter(const double& p);

    ///Set/Get the SeparationAlgorithm HistogramBins parameter
    const int GetSeparationAlgorithmHistogramBinsParameter() const;
    ///Set/Get the SeparationAlgorithm HistogramBins parameter
    bool SetSeparationAlgorithmHistogramBinsParameter(const int& p);

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
    //The analysis model tag and one attribute
    static inline const char* analysisModelTag() { return "analysis-model"; }
    static inline const char* analysisModelNameAttribute() { return "model-name"; }
    //The algorithm tag and one attribute
    static inline const char* algorithmTag() { return "algorithm"; }
    static inline const char* algorithmNameAttribute() { return "alg-name"; }
    //The parameter tag and one attribute
    static inline const char* parameterTag() { return "parameter"; }
    static inline const char* parameterTypeAttribute() { return "param-type"; }

    //Parameter type possible values
    static inline const char* pTypeNumPixels() { return "num-pixels"; }
    static inline const char* pTypeThreshold() { return "threshold"; }
    static inline const char* pTypePercentile() { return "percentile"; }
    static inline const char* pTypeHistoBins() { return "histo-bins"; }

private:
    ///Build the XMLDocument data structure
    bool BuildXMLDocument();
    ///Check if the basic structure of the member XMLDocument has been assembled
    bool CheckXMLDocument();
    ///Clear the entire contents of the member XMLDocument
    bool ClearXMLDocument();

    ///If file is able to be opened, write the current stain profile to file as XML
    tinyxml2::XMLError writeStainProfileToXML(const std::string &);
    ///If file is able to be opened, read from an XML file and fill variables in this class
    tinyxml2::XMLError readStainProfileFromXMLFile(const std::string &);
    ///Read a stain profile from a string
    tinyxml2::XMLError readStainProfileFromXMLString(const char *, size_t);
    //Parse XML document
    tinyxml2::XMLError parseXMLDoc();

private:
    ///Store the list of possible stain analysis model names here
    const std::vector<std::string> m_stainAnalysisModelOptions;
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
    ///stain analysis model element
    tinyxml2::XMLElement* m_analysisModelElement;
    ///algorithm element
    tinyxml2::XMLElement* m_algorithmElement;
};

#endif
