# Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
#
# File setting up CPack-generator specific options.
#

if( "${CPACK_GENERATOR}" STREQUAL "PackageMaker" )
   set( CPACK_PACKAGING_INSTALL_PREFIX "/${CPACK_PACKAGE_INSTALL_DIRECTORY}" )
   set( CPACK_PACKAGE_DEFAULT_LOCATION "/Applications/ATLAS" )
elseif( "${CPACK_GENERATOR}" STREQUAL "TGZ" )
   set( CPACK_INCLUDE_TOPLEVEL_DIRECTORY FALSE )
   set( CPACK_PACKAGING_INSTALL_PREFIX "/${CPACK_PACKAGE_INSTALL_DIRECTORY}" )
   set( CPACK_SET_DESTDIR TRUE )
endif()
