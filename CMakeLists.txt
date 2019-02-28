PROJECT( StainAnalysis )
CMAKE_MINIMUM_REQUIRED( VERSION 3.8 )
## Set the C++ standard to use
SET(CMAKE_CXX_STANDARD 17)

##
## Load the Sedeen dependencies
SET(PROGRAMFILESX86 "PROGRAMFILES\(X86\)")
FIND_PACKAGE( SEDEENSDK REQUIRED 
                HINTS ../../.. 
				   "$ENV{PROGRAMFILESX86}/Sedeen Viewer SDK/v5.4.0.20181217/msvc2017"
                   "$ENV{PROGRAMFILES}/Sedeen Viewer SDK/v5.4.0.20181217/msvc2017" )

##
## The ColorDeconvolutionKernel depends on boost.filesystem
IF(NOT BOOST_ROOT)
  SET(BOOST_ROOT "BOOST_ROOT-NOTFOUND" CACHE PATH "Preferred installation prefix of the Boost C++ library")
ENDIF()
IF(NOT BOOST_VERSION)
  SET(BOOST_VERSION "BOOST_VERSION-NOTFOUND" CACHE STRING "Boost library version number")
ENDIF()
FIND_PACKAGE(Boost ${BOOST_VERSION} REQUIRED COMPONENTS
               filesystem )

INCLUDE_DIRECTORIES( ${INCLUDE_DIRECTORIES} "${SEDEENSDK_INCLUDE_DIR}" "${BOOST_ROOT}")

LINK_DIRECTORIES( ${LINK_DIRECTORIES} "${SEDEENSDK_LIBRARY_DIR}" "${BOOST_ROOT}/stage/lib")

##
## Build the code into a module library
ADD_LIBRARY( StainAnalysis MODULE StainAnalysis.cpp StainAnalysis.h ColorDeconvolutionKernel.cpp ColorDeconvolutionKernel.h)

# Link the library against the Sedeen libraries
# NOTE: The QT libraries must be linked first.
TARGET_LINK_LIBRARIES( StainAnalysis ${SEDEENSDK_LIBRARIES} 
					     Boost::filesystem )

##
## Install the plugin in the sedeen plugins directory
IF( ${PATHCORE_FOUND} )
  INSTALL( TARGETS StainAnalysis 
           LIBRARY DESTINATION "${PATHCORE_DIR}/plugins" )
ENDIF()
