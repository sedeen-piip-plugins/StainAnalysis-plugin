PROJECT( StainAnalysis-plugin )
CMAKE_MINIMUM_REQUIRED( VERSION 3.8 )
# Enable C++17 features
SET(CMAKE_CXX_STANDARD 17)

# Define project description variables
SET( DISPLAY_NAME_TEXT "${PROJECT_NAME}" CACHE STRING "Name of the plugin as it should be displayed in Sedeen Viewer")
SET( SUPPORT_URL_TEXT "http://pathcore.com/support/plugin/info/${PROJECT_NAME}" CACHE STRING "Location users can find help with the plugin" )
SET( DEVELOPER_TEXT "Sunnybrook Research Institute" CACHE STRING "Name of the author or organization that created the plugin" )

# Load the Sedeen dependencies
SET(PROGRAMFILESX86 "PROGRAMFILES\(X86\)")
FIND_PACKAGE( SEDEENSDK REQUIRED 
                HINTS ../../.. 
                  "$ENV{${PROGRAMFILESX86}}/Sedeen Viewer SDK/v5.4.0.20181217/msvc2017"
                  "$ENV{PROGRAMFILES}/Sedeen Viewer SDK/v5.4.0.20181217/msvc2017" )

# The ColorDeconvolutionKernel depends on boost.filesystem
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

# Build the code into a module library
ADD_LIBRARY( ${PROJECT_NAME} MODULE ${PROJECT_NAME}.cpp ${PROJECT_NAME}.h ColorDeconvolutionKernel.cpp ColorDeconvolutionKernel.h)

# Link the library against the Sedeen SDK libraries
TARGET_LINK_LIBRARIES( ${PROJECT_NAME} ${SEDEENSDK_LIBRARIES}
					     Boost::filesystem )

#Create or update the .info file in the build directory
STRING( TIMESTAMP DATE_CREATED_TEXT "%Y-%m-%d" )
CONFIGURE_FILE( "${PROJECT_NAME}.info.in" "${PROJECT_NAME}.info" )

# Set the install destination directory
IF( NOT PLUGIN_DESTINATION_DIR )
  IF( ${SEDEEN_FOUND} )
	SET( TEMPPLUGINDIR "${PATHCORE_DIR}/plugins/cpp/piip/${PROJECT_NAME}" )
  ELSE()
	SET( TEMPPLUGINDIR "PLUGIN_DESTINATION_DIR-NOTFOUND" )
	MESSAGE( SEND_ERROR "PLUGIN_DESTINATION_DIR not found. Set this to the target installation directory of the plugin within Sedeen Viewer (e.g. $ENV{PROGRAMFILES}/Sedeen Viewer/plugins/cpp/piip/${PROJECT_NAME}).")
  ENDIF()
  SET(PLUGIN_DESTINATION_DIR ${TEMPPLUGINDIR} CACHE PATH "Installation directory for the plugin within Sedeen Viewer")
ENDIF()

# Install the plugin and .info file in the PLUGIN_DESTINATION_DIR directory
IF( ${SEDEEN_FOUND} )
  install(TARGETS ${PROJECT_NAME}
          LIBRARY DESTINATION "${PLUGIN_DESTINATION_DIR}")
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.info"
          DESTINATION "${PLUGIN_DESTINATION_DIR}")
ENDIF()