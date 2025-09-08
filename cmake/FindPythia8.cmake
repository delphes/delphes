# Find the Pythia8 includes and library.
#
# This module defines
# PYTHIA8_INCLUDE_DIR   where to locate Pythia.h file
# PYTHIA8_LIBRARY       where to find the libpythia8 library
# PYTHIA8_LIBRARIES     (not cached) the libraries to link against to use Pythia8
# PYTHIA8_FOUND         if false, you cannot build anything that requires Pythia8
# PYTHIA8_VERSION       version of Pythia8 if found

set(_pythia8dirs
    ${PYTHIA8}
    $ENV{PYTHIA8}
    ${PYTHIA8_DIR}
    $ENV{PYTHIA8_DIR}
    /usr
    /usr/local
    /opt/pythia
    /opt/pythia8)

find_path(PYTHIA8_INCLUDE_DIR
          NAMES Pythia8/Pythia.h
          HINTS ${_pythia8dirs}
          PATH_SUFFIXES include include/Pythia8 include/pythia8)

find_library(PYTHIA8_LIBRARY
             NAMES pythia8 Pythia8
             HINTS ${_pythia8dirs}
             PATH_SUFFIXES lib)

unset(_pythia8dirs)

set(PYTHIA8_INCLUDE_DIRS ${PYTHIA8_INCLUDE_DIR} ${PYTHIA8_INCLUDE_DIR}/Pythia8)
set(PYTHIA8_LIBRARIES ${PYTHIA8_LIBRARY})

# handle the QUIETLY and REQUIRED arguments and set PYTHIA8_FOUND to TRUE if
# all listed variables are TRUE

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Pythia8 DEFAULT_MSG PYTHIA8_INCLUDE_DIR PYTHIA8_LIBRARY)
mark_as_advanced(PYTHIA8_INCLUDE_DIR PYTHIA8_LIBRARY)
