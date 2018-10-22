# Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
#
# This file describes to CTest how to build an entire release, and upload
# the build results to a CDash server.
#

# Set some basic properties for the build:
set( CTEST_SOURCE_DIRECTORY /global/u1/t/tongou/WSInspector/source )
set( CTEST_BINARY_DIRECTORY /global/u1/t/tongou/WSInspector/build )

set( CTEST_CUSTOM_MAXIMUM_NUMBER_OF_ERRORS   1000000 )
set( CTEST_CUSTOM_MAXIMUM_NUMBER_OF_WARNINGS 1000000 )

set( CTEST_COMMAND "/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/Cmake/3.11.0/Linux-x86_64/bin/ctest -D Nightly" )
set( CTEST_CMAKE_COMMAND /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/Cmake/3.11.0/Linux-x86_64/bin/cmake )

# The nightly build always uses launchers to capture the full build logs
# of the packages.
set( CTEST_USE_LAUNCHERS 1 )

# Set how the build should be done exactly:
set( CTEST_CMAKE_GENERATOR "Unix Makefiles" )
set( CTEST_BUILD_CONFIGURATION RelWithDebInfo )
if( ${CTEST_CMAKE_GENERATOR} STREQUAL "Unix Makefiles" )
   set( CTEST_BUILD_FLAGS -j32 )
endif()

# Identify the build:
set( CTEST_BUILD_NAME x86_64-slc6-gcc62-opt )
set( CTEST_SITE pdsf7 )

# Start a nightly build:
message( "Starting nightly build for project WorkDir..." )
ctest_start( Nightly )

# Find out the tag used for the current build:
file( READ ${CTEST_BINARY_DIRECTORY}/Testing/TAG _tagFile )
set( _tag )
if( _tagFile MATCHES "([0-9]*-[0-9]*).*" )
   set( _tag ${CMAKE_MATCH_1} )
endif()
message( "Using nightly tag: ${_tag}" )

# Find the log extractor scripts:
find_file( _buildExtractor NAMES extractCTestBuild.py
   PATH_SUFFIXES scripts PATHS /cvmfs/atlas.cern.ch/repo/sw/software/21.2/AnalysisBase/21.2.27/InstallArea/x86_64-slc6-gcc62-opt/cmake/modules;/cvmfs/atlas.cern.ch/repo/sw/software/21.2/AnalysisBaseExternals/21.2.27/InstallArea/x86_64-slc6-gcc62-opt/cmake/modules )
if( _buildExtractor MATCHES "NOTFOUND" )
   message( FATAL_ERROR "Couldn't find the extractCTestBuild.py script" )
endif()
message( "Using: ${_buildExtractor}" )
find_file( _testExtractor NAMES extractCTestTest.py
   PATH_SUFFIXES scripts PATHS /cvmfs/atlas.cern.ch/repo/sw/software/21.2/AnalysisBase/21.2.27/InstallArea/x86_64-slc6-gcc62-opt/cmake/modules;/cvmfs/atlas.cern.ch/repo/sw/software/21.2/AnalysisBaseExternals/21.2.27/InstallArea/x86_64-slc6-gcc62-opt/cmake/modules )
if( _testExtractor MATCHES "NOTFOUND" )
   message( FATAL_ERROR "Couldn't find the extractCTestTest.py script" )
endif()
message( "Using: ${_testExtractor}" )

# Upload the project description file to CDash:
if( EXISTS ${CTEST_BINARY_DIRECTORY}/Project.xml )
   message( "Uploading the sub-project description file" )
   ctest_submit( FILES ${CTEST_BINARY_DIRECTORY}/Project.xml )
endif()

# Execute the project configuration just once:
message( "Configuring the build..." )
ctest_configure( BUILD ${CTEST_BINARY_DIRECTORY}
   SOURCE ${CTEST_SOURCE_DIRECTORY}
   OPTIONS "-DCTEST_USE_LAUNCHERS=${CTEST_USE_LAUNCHERS}" )
message( "Submitting the configuration results" )
ctest_submit( PARTS Configure )

# The subprojects to build/test:
set( _subprojects  )
list( LENGTH _subprojects _nSubprojects )

# Run the full build in one go:
message( "Building the entire project (${_nSubprojects} packages)..." )
ctest_build( BUILD ${CTEST_BINARY_DIRECTORY} APPEND )

# Make a copy of the build file produced:
message( "Making a copy of the full build log..." )
execute_process( COMMAND ${CTEST_CMAKE_COMMAND} -E copy
   ${CTEST_BINARY_DIRECTORY}/Testing/${_tag}/Build.xml
   ${CTEST_BINARY_DIRECTORY}/Testing/${_tag}/FullBuild.xml )

# If there are no subprojects, upload the build results now, with the
# full build log(s):
if( NOT _subprojects )
   set( _toSubmit Build )
   file( GLOB _buildLogs ${CTEST_BINARY_DIRECTORY}/BuildLogs/*.log )
   if( _buildLogs )
      set( CTEST_NOTES_FILES ${_buildLogs} )
      list( APPEND _toSubmit Notes )
   endif()
   ctest_submit( PARTS ${_toSubmit} )
endif()

# Run all the tests in one go:
message( "Running the tests of all the packages..." )
ctest_test( BUILD ${CTEST_BINARY_DIRECTORY}
   PARALLEL_LEVEL 32 )

# Make a copy of the test file produced:
message( "Making a copy of the full test log..." )
execute_process( COMMAND ${CTEST_CMAKE_COMMAND} -E copy
   ${CTEST_BINARY_DIRECTORY}/Testing/${_tag}/Test.xml
   ${CTEST_BINARY_DIRECTORY}/Testing/${_tag}/FullTest.xml )

# If there are no subprojects, upload the test results now:
if( NOT _subprojects )
   ctest_submit( PARTS Test )
endif()

# Loop over all defined subprojects, and upload their build and test results:
foreach( _subproject ${_subprojects} )

   # Let the user know what's happening:
   message( "Uploading results for package: ${_subproject}" )

   # Define which subproject we're looking at:
   set_property( GLOBAL PROPERTY SubProject ${_subproject} )
   set_property( GLOBAL PROPERTY Label ${_subproject} )

   # The list of components to submit:
   set( _toSubmit Build )

   # The log file that was created for the build, if it exists:
   if( EXISTS ${CTEST_BINARY_DIRECTORY}/BuildLogs/${_subproject}.log )
      set( CTEST_NOTES_FILES
         ${CTEST_BINARY_DIRECTORY}/BuildLogs/${_subproject}.log )
      list( APPEND _toSubmit Notes )
   endif()

   # Extract the build results:
   execute_process( COMMAND ${_buildExtractor}
      -i ${CTEST_BINARY_DIRECTORY}/Testing/${_tag}/FullBuild.xml
      -o ${CTEST_BINARY_DIRECTORY}/Testing/${_tag}/Build.xml
      -l ${_subproject} )
   # And also save them for debugging:
   execute_process( COMMAND ${CTEST_CMAKE_COMMAND} -E copy
      ${CTEST_BINARY_DIRECTORY}/Testing/${_tag}/Build.xml
      ${CTEST_BINARY_DIRECTORY}/Testing/${_tag}/${_subproject}Build.xml )

   # Submit the build results:
   ctest_submit( PARTS ${_toSubmit} )

   # Extract the test results:
   execute_process( COMMAND ${_testExtractor}
      -i ${CTEST_BINARY_DIRECTORY}/Testing/${_tag}/FullTest.xml
      -o ${CTEST_BINARY_DIRECTORY}/Testing/${_tag}/Test.xml
      -l ${_subproject} )
   # And also save them for debugging:
   execute_process( COMMAND ${CTEST_CMAKE_COMMAND} -E copy
      ${CTEST_BINARY_DIRECTORY}/Testing/${_tag}/Test.xml
      ${CTEST_BINARY_DIRECTORY}/Testing/${_tag}/${_subproject}Test.xml )

   # Submit the build results:
   ctest_submit( PARTS Test )

endforeach()
