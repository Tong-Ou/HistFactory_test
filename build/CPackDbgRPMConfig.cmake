# Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
#
# Configuration for building an RPM file out of the .dbg files produced during
# a "RelWithDebInfo" build.
#

# Set some basic variables:
set( CPACK_CMAKE_GENERATOR "Unix Makefiles" )
set( CPACK_GENERATOR "RPM" )
set( CPACK_INSTALL_CMAKE_PROJECTS "" )
set( CPACK_PACKAGE_DESCRIPTION_FILE "/global/u1/t/tongou/WSInspector/build/README.txt" )
set( CPACK_RESOURCE_FILE_README "/global/u1/t/tongou/WSInspector/build/README.txt" )
set( CPACK_SET_DESTDIR "OFF" )
set( CPACK_OUTPUT_CONFIG_FILE "/global/u1/t/tongou/WSInspector/build/CPackConfig.cmake" )

# Set the package name, version, and other minor information:
set( CPACK_PACKAGE_NAME "WorkDir" )
set( CPACK_PACKAGE_VERSION "1.0.0" )
set( CPACK_PACKAGE_VERSION_MAJOR "1" )
set( CPACK_PACKAGE_VERSION_MINOR "0" )
set( CPACK_PACKAGE_VERSION_PATCH "0" )
set( CPACK_PACKAGE_VENDOR "ATLAS Collaboration" )
set( CPACK_PACKAGE_CONTACT "atlas-sw-core@cern.ch" )

# Set the RPM file name:
set( CPACK_PACKAGE_FILE_NAME
   "WorkDir_${CPACK_PACKAGE_VERSION}_x86_64-slc6-gcc62-opt-dbg" )

# Set up where the project is installed:
set( CPACK_PACKAGE_RELOCATABLE TRUE )
set( CPACK_PACKAGE_DEFAULT_LOCATION "/usr" )
set( CPACK_PACKAGE_INSTALL_DIRECTORY
   "WorkDir/${CPACK_PACKAGE_VERSION}/InstallArea/x86_64-slc6-gcc62-opt" )
set( CPACK_INSTALL_PREFIX
   "${CPACK_PACKAGE_DEFAULT_LOCATION}/${CPACK_PACKAGE_INSTALL_DIRECTORY}" )

# Set up how the debug files should be installed for the RPM build:
set( CPACK_INSTALL_COMMANDS
   "/global/u1/t/tongou/WSInspector/build/CMakeFiles/cpack_dbg_install.sh" )

# Set the RPM specific variables:
set( CPACK_RPM_PACKAGE_NAME "${CPACK_PACKAGE_FILE_NAME}" )
set( CPACK_RPM_PACKAGE_VERSION "${CPACK_PACKAGE_VERSION}" )
set( CPACK_RPM_PACKAGE_PROVIDES "/bin/sh" )
set( CPACK_RPM_PACKAGE_ARCHITECTURE "noarch" )
set( CPACK_RPM_PACKAGE_GROUP "ATLAS Software" )
set( CPACK_RPM_PACKAGE_LICENSE "Simplified BSD" )
set( CPACK_RPM_SPEC_MORE_DEFINE "%global __os_install_post %{nil}" )

# Make the RPM depend on the main (opt) RPM:
set( CPACK_RPM_PACKAGE_AUTOREQ " no" )
set( CPACK_RPM_PACKAGE_AUTOREQPROV " no" )
set( CPACK_RPM_PACKAGE_REQUIRES
   "WorkDir_${CPACK_PACKAGE_VERSION}_x86_64-slc6-gcc62-opt" )
