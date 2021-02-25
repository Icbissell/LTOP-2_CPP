# This file builds some boost components using the boost directory in this repo
include( ExternalProject )

set( boost_INSTALL ${CMAKE_CURRENT_BINARY_DIR}/third_party/boost )
set( boost_INCLUDE_DIR ${boost_INSTALL}/include )
set( boost_LIB_DIR ${boost_INSTALL}/lib )

file( MAKE_DIRECTORY ${boost_INCLUDE_DIR} )

ExternalProject_Add( external_boost
        PREFIX boost
		  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/LTOP2-CPP/boost_1_73_0
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND find . -name "*.sh" | xargs chmod a+x && ./bootstrap.sh
            --with-libraries=filesystem
            --with-libraries=iostreams
            --with-libraries=system
            --prefix=<INSTALL_DIR>
        BUILD_COMMAND
        ./b2 install link=static variant=release threading=multi runtime-link=static
        INSTALL_COMMAND ""
        INSTALL_DIR ${boost_INSTALL} )

set( boost_LIBRARY_SUFFIX .a )

add_library( boost::iostreams STATIC IMPORTED )
set_property( TARGET boost::iostreams PROPERTY IMPORTED_LOCATION ${boost_LIB_DIR}/libboost_iostreams${boost_LIBRARY_SUFFIX} )
set_property( TARGET boost::iostreams PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${boost_INCLUDE_DIR} )
add_dependencies( boost::iostreams external_boost )

add_library( boost::system STATIC IMPORTED )
set_property( TARGET boost::system PROPERTY IMPORTED_LOCATION ${boost_LIB_DIR}/libboost_system${boost_LIBRARY_SUFFIX} )
set_property( TARGET boost::system PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${boost_INCLUDE_DIR} )
add_dependencies( boost::system external_boost )

add_library( boost::filesystem STATIC IMPORTED )
set_property( TARGET boost::filesystem PROPERTY IMPORTED_LOCATION ${boost_LIB_DIR}/libboost_filesystem${boost_LIBRARY_SUFFIX} )
set_property( TARGET boost::filesystem PROPERTY INTERFACE_LINK_LIBRARIES boost::system )
set_property( TARGET boost::filesystem PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${boost_INCLUDE_DIR} )
add_dependencies( boost::filesystem external_boost )
