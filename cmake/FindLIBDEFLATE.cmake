if(UNIX)
    find_package(PkgConfig QUIET)
    # pkg-config support added in libdeflate v1.9
    pkg_check_modules(_LIBDEFLATE QUIET libdeflate)
else()
    message(FATAL_ERROR "BISCUIT only works on Unix systems")
endif()

find_path(LIBDEFLATE_INCLUDE_DIR
    NAMES libdeflate.h
    HINTS ${_LIBDEFLATE_INCLUDEDIR})
find_library(LIBDEFLATE_LIBRARY
    NAMES deflate
    HINTS ${_LIBDEFLATE_LIBDIR})

set(LIBDEFLATE_INCLUDE_DIRS ${LIBDEFLATE_INCLUDE_DIR})
set(LIBDEFLATE_LIBRARIES ${LIBDEFLATE_LIBRARY})

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
