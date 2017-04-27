PROJECT( StainAnalysis )
CMAKE_MINIMUM_REQUIRED( VERSION 2.8 )
CMAKE_POLICY( SET CMP0053 OLD )

##
## Load the Sedeen dependencies
FIND_PACKAGE( SEDEENSDK REQUIRED 
                HINTS ../../.. 
                "$ENV{PROGRAMFILES(X86)}/Sedeen Viewer SDK/v5.2.1.427/msvc2012"
                "$C/Azadeh/Sedeen Viewer SDK/v5.2.1.427/msvc2012" )


INCLUDE_DIRECTORIES( "${SEDEENSDK_INCLUDE_DIR}" )

LINK_DIRECTORIES( "${SEDEENSDK_LIBRARY_DIR}" )

##
## Build the code into a module library
ADD_LIBRARY( StainAnalysis MODULE StainAnalysis.cpp StainAnalysis.h ColorDeconvolutionKernel.cpp ColorDeconvolutionKernel.h)

# Link the library against the Sedeen libraries
# NOTE: The QT libraries must be linked first.
TARGET_LINK_LIBRARIES( StainAnalysis ${PATHCORE_LIBRARIES} )

##
## Install the plugin in the sedeen plugins directory
IF( ${PATHCORE_FOUND} )
  INSTALL( TARGETS StainAnalysis 
           LIBRARY DESTINATION "${PATHCORE_DIR}/plugins" )
ENDIF()
