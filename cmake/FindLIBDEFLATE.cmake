# Search for libdeflate using pkg-config
if(UNIX)
    find_package(PkgConfig QUIET)
    pkg_check_modules(_LIBDEFLATE QUIET libdeflate)
else()
    message(FATAL_ERROR "BISCUIT currently only works on Unix systems")
endif()

# Get include directory and shared library path
find_path(LIBDEFLATE_INCLUDE_DIR
    NAMES libdeflate.h
    HINTS ${_LIBDEFLATE_INCLUDEDIR}
    )
find_library(LIBDEFLATE_LIBRARY
    NAMES deflate
    HINTS ${_LIBDEFLATE_LIBDIR}
    )

# Set variables for use in CMakeLists files
set(LIBDEFLATE_INCLUDE_DIRS ${LIBDEFLATE_INCLUDE_DIR})
set(LIBDEFLATE_LIBRARIES ${LIBDEFLATE_LIBRARY})

# Pull version from pkg-config if possible, otherwise pull from libdeflate header
if(_LIBDEFLATE_VERSION)
    set(LIBDEFLATE_VERSION ${_LIBDEFLATE_VERSION})
elseif(LIBDEFLATE_INCLUDE_DIR)
    file(STRINGS "${LIBDEFLATE_INCLUDE_DIR}/libdeflate.h" LIBDEFLATE_VERSION_STR
        REGEX "^#define[\t ]+LIBDEFLATE_VERSION_STRING[\t ]+\"[^\"]+\"")
    if(LIBDEFLATE_VERSION_STR MATCHES "\"([^\"]+)\"")
        set(LIBDEFLATE_VERSION "${CMAKE_MATCH_1}")
    endif()
endif()

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(LIBDEFLATE
    REQUIRED_VARS
        LIBDEFLATE_INCLUDE_DIR
        LIBDEFLATE_LIBRARY
    VERSION_VAR LIBDEFLATE_VERSION)

mark_as_advanced(LIBDEFLATE_INCLUDE_DIR LIBDEFLATE_LIBRARY)
