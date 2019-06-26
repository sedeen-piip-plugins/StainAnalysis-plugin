PROJECT( StainAnalysis-plugin )
CMAKE_MINIMUM_REQUIRED( VERSION 3.8 )
# Enable C++17 features
SET(CMAKE_CXX_STANDARD 17)

# Define project description variables
SET( DISPLAY_NAME_TEXT "Stain Analysis 5.4.1.20190327" CACHE STRING "Name of the plugin as it should be displayed in Sedeen Viewer")
SET( SUPPORT_URL_TEXT "http://pathcore.com/support/plugin/info/${PROJECT_NAME}" CACHE STRING "Location users can find help with the plugin" )
SET( DEVELOPER_TEXT "Sunnybrook Research Institute" CACHE STRING "Name of the author or organization that created the plugin" )

# Load the Sedeen dependencies
SET(PROGRAMFILESX86 "PROGRAMFILES\(X86\)")
FIND_PACKAGE( SEDEENSDK REQUIRED 
                HINTS ../../.. 
                  "$ENV{${PROGRAMFILESX86}}/Sedeen Viewer SDK/v5.4.1.20190327/msvc2017"
                  "$ENV{PROGRAMFILES}/Sedeen Viewer SDK/v5.4.1.20190327/msvc2017" )

#Add OpenMP for parallel programming
#This is not essential to the plugin and can be removed
FIND_PACKAGE(OPENMP)
IF(OPENMP_FOUND)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
ENDIF()

IF(NOT BOOST_ROOT)
  SET(BOOST_ROOT "BOOST_ROOT-NOTFOUND" CACHE PATH "Preferred installation prefix of the Boost C++ library")
ENDIF()
IF(NOT BOOST_VERSION)
  SET(BOOST_VERSION "BOOST_VERSION-NOTFOUND" CACHE STRING "Boost library version number")
ENDIF()
FIND_PACKAGE(Boost ${BOOST_VERSION} )

#Inlude TinyXML2 for XML file reading and writing
IF(NOT TinyXML2_INCLUDE_DIR)
  SET(TinyXML2_INCLUDE_DIR "TinyXML2_INCLUDE_DIR-NOTFOUND" CACHE PATH "Location of include files for TinyXML2")
ENDIF()
IF(NOT TinyXML2_LIB_DIR)
  SET(TinyXML2_LIB_DIR "TinyXML2_LIB_DIR-NOTFOUND" CACHE PATH "Location of library files for TinyXML2")
ENDIF()
FIND_PACKAGE(TinyXML2)

INCLUDE_DIRECTORIES( ${INCLUDE_DIRECTORIES} ${SEDEENSDK_INCLUDE_DIR} ${BOOST_ROOT} ${TinyXML2_INCLUDE_DIR})
LINK_DIRECTORIES( ${LINK_DIRECTORIES} ${SEDEENSDK_LIBRARY_DIR} ${TinyXML2_LIB_DIR})

# Build the code into a module library
ADD_LIBRARY( ${PROJECT_NAME} MODULE 
                 ${PROJECT_NAME}.cpp ${PROJECT_NAME}.h 
                 StainProfile.h StainProfile.cpp 
                 ColorDeconvolutionKernel.h ColorDeconvolutionKernel.cpp
				 StainVectorMath.h StainVectorMath.cpp
                 )

# Link the library against the Sedeen SDK libraries
TARGET_LINK_LIBRARIES( ${PROJECT_NAME} ${SEDEENSDK_LIBRARIES}
                         TinyXML2 )

#Create or update the .info file in the build directory
STRING( TIMESTAMP DATE_CREATED_TEXT "%Y-%m-%d" )
CONFIGURE_FILE( "infoTemplate.info.in" "${PROJECT_NAME}.info" )

#Set the relative directory where the plugin should be located
SET( PLUGIN_RELATIVE_DIR "plugins/cpp/piip/${PROJECT_NAME}" )
#Add to preprocessor definitions
ADD_DEFINITIONS( -DPLUGIN_RELATIVE_DIR=${PLUGIN_RELATIVE_DIR} )

# Set the install destination directory
IF( NOT PLUGIN_DESTINATION_DIR )
  IF( ${SEDEEN_FOUND} )
	SET( TEMPPLUGINDIR "${PATHCORE_DIR}/${PLUGIN_RELATIVE_DIR}" )
  ELSE()
	SET( TEMPPLUGINDIR "PLUGIN_DESTINATION_DIR-NOTFOUND" )
	MESSAGE( SEND_ERROR "PLUGIN_DESTINATION_DIR not found. Set this to the target installation directory of the plugin within Sedeen Viewer (e.g. $ENV{PROGRAMFILES}/Sedeen Viewer/plugins/cpp/piip/${PROJECT_NAME}).")
  ENDIF()
  SET(PLUGIN_DESTINATION_DIR ${TEMPPLUGINDIR} CACHE PATH "Installation directory for the plugin within Sedeen Viewer")
ENDIF()

# Install the plugin and .info file in the PLUGIN_DESTINATION_DIR directory
IF( ${SEDEEN_FOUND} )
  INSTALL(TARGETS ${PROJECT_NAME}
          LIBRARY DESTINATION "${PLUGIN_DESTINATION_DIR}")
  INSTALL(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.info"
          DESTINATION "${PLUGIN_DESTINATION_DIR}")
  INSTALL(FILES "defaultprofiles/HematoxylinPEosin.xml"
                "defaultprofiles/HematoxylinPDAB.xml"
                "defaultprofiles/HematoxylinPEosinPDAB.xml"
            DESTINATION "${PLUGIN_DESTINATION_DIR}")
#  INSTALL(FILES 
#            DESTINATION "${PLUGIN_DESTINATION_DIR}")

ENDIF()