PROJECT( StainAnalysis )
CMAKE_MINIMUM_REQUIRED( VERSION 2.8 )
CMAKE_POLICY( SET CMP0053 OLD )

##
## Load the Sedeen dependencies
FIND_PACKAGE( SEDEENSDK REQUIRED 
                HINTS ../../.. 
                "$ENV{PROGRAMFILES(X86)}/Sedeen Viewer SDK/v5.2.1.427/msvc2012"
                "$C/Azadeh/Sedeen Viewer SDK/v5.2.1.427/msvc2012" )

##
## The ColorDeconvolutionKernel depends on boost.filesystem
set(BOOST_ROOT "" CACHE PATH "Preferred installation prefix of the Boost C++ library")
#set(BOOST_INCLUDEDIR "${BOOST_ROOT}/include" CACHE PATH "Preferred include directory of the Boost C++ library e.g. <prefix>/include")
#set(BOOST_LIBRARYDIR "${BOOST_ROOT}/libs" CACHE PATH "Preferred library directory of the Boost C++ library e.g. <prefix>/lib")

#set(Boost_LIBRARY_DIR_DEBUG "" CACHE PATH "Directory containing debug Boost libraries")
#set(Boost_LIBRARY_DIR_RELEASE "" CACHE PATH "Directory containing release Boost libraries")

#"The directory containing a CMake configuration file for Boost"
#"Boost filesystem directory (debug)"
#"Boost filesystem directory (release)"

FIND_PACKAGE(Boost 1.69 REQUIRED COMPONENTS
               filesystem )

INCLUDE_DIRECTORIES( "${SEDEENSDK_INCLUDE_DIR}" )

LINK_DIRECTORIES( "${SEDEENSDK_LIBRARY_DIR}" )

##
## Build the code into a module library
ADD_LIBRARY( StainAnalysis MODULE StainAnalysis.cpp StainAnalysis.h ColorDeconvolutionKernel.cpp ColorDeconvolutionKernel.h)

# Link the library against the Sedeen libraries
# NOTE: The QT libraries must be linked first.
TARGET_LINK_LIBRARIES( StainAnalysis ${PATHCORE_LIBRARIES} 
					     Boost::filesystem )

##
## Install the plugin in the sedeen plugins directory
IF( ${PATHCORE_FOUND} )
  INSTALL( TARGETS StainAnalysis 
           LIBRARY DESTINATION "${PATHCORE_DIR}/plugins" )
ENDIF()
