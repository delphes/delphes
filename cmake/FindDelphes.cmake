# Find the Delphes includes and library.

set(_delphesdirs
    ${DELPHES}
    $ENV{DELPHES}
    ${DELPHES_DIR}
    $ENV{DELPHES_DIR}
    /usr
    /usr/local
    /opt/delphes)

find_library(DELPHES_LIBRARY
             NAMES Delphes delphes
             HINTS ${_delphesdirs}
             PATH_SUFFIXES lib)

find_path(DELPHES_INCLUDE_DIR
          NAMES classes/DelphesClasses.h modules/Delphes.h
          HINTS ${_delphesdirs}
          PATH_SUFFIXES include)

find_path(DELPHES_EXTERNALS_INCLUDE_DIR
          NAMES ExRootAnalysis/ExRootConfReader.h ExRootAnalysis/ExRootTreeReader.h
          HINTS ${_delphesdirs}
          PATH_SUFFIXES include)

find_path(DELPHES_BINARY_DIR
          NAMES DelphesROOT DelphesSTDHEP
          HINTS ${_delphesdirs}
          PATH_SUFFIXES bin)

find_path(DELPHES_CARDS_DIR
          NAMES delphes_card_ATLAS.tcl delphes_card_CMS.tcl
          HINTS ${_delphesdirs}
          PATH_SUFFIXES cards)

unset(_delphesdirs)

set(DELPHES_INCLUDE_DIRS ${DELPHES_INCLUDE_DIR} ${DELPHES_EXTERNALS_INCLUDE_DIR})
set(DELPHES_EXTERNALS_INCLUDE_DIRS ${DELPHES_EXTERNALS_INCLUDE_DIR})
set(DELPHES_LIBRARIES ${DELPHES_LIBRARY})

# Delphes does not offer an obvious version indicator, but we need to know
# whether the TrackCovariance module is available or not. So here we simply
# check whether the corresponding header is installed
find_file(DELPHES_TRACK_COV_HEADER modules/TrackCovariance.h PATHS ${DELPHES_INCLUDE_DIRS} NO_DEFAULT_PATHS)

# handle the QUIETLY and REQUIRED arguments and set DELPHES_FOUND to TRUE
# if all listed variables are TRUE

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Delphes DEFAULT_MSG DELPHES_INCLUDE_DIR DELPHES_EXTERNALS_INCLUDE_DIR DELPHES_LIBRARY)
mark_as_advanced(DELPHES_INCLUDE_DIR DELPHES_EXTERNALS_INCLUDE_DIR DELPHES_LIBRARY DELPHES_BINARY_DIR DELPHES_TRACK_COV_HEADER)
