# Find an installed FastJet contrib with the contribs used by Delphes
# (Nsubjettiness, ValenciaPlugin, RecursiveTools)
#
# Hints: FJCONTRIB_ROOT (CMake or environment variable), the FastJet location
#
# Uses the shared libfastjetcontribfragile if available (make fragile-shared-install),
# otherwise the individual contrib libraries (which then need to be built with -fPIC).
#
# Defines
#   FastJetContrib_FOUND, FastJetContrib_INCLUDE_DIRS, FastJetContrib_LIBRARIES
#   the imported target FastJet::Contrib

set(_fjcontrib_hints ${FJCONTRIB_ROOT} $ENV{FJCONTRIB_ROOT})
if(FastJet_INCLUDE_DIR)
  get_filename_component(_fjcontrib_fastjet_prefix ${FastJet_INCLUDE_DIR} DIRECTORY)
  list(APPEND _fjcontrib_hints ${_fjcontrib_fastjet_prefix})
endif()

set(_fjcontrib_headers Nsubjettiness.hh NjettinessPlugin.hh ExtraRecombiners.hh ValenciaPlugin.hh SoftDrop.hh)

find_path(_fjcontrib_header_dir fastjet/contrib/Nsubjettiness.hh
  HINTS ${_fjcontrib_hints} PATH_SUFFIXES include)

if(_fjcontrib_header_dir)
  get_filename_component(_fjcontrib_header ${_fjcontrib_header_dir}/fastjet/contrib/Nsubjettiness.hh REALPATH)
  get_filename_component(_fjcontrib_header_dir ${_fjcontrib_header} DIRECTORY)
  get_filename_component(_fjcontrib_header_dir ${_fjcontrib_header_dir} DIRECTORY)
  get_filename_component(FastJetContrib_INCLUDE_DIR ${_fjcontrib_header_dir} DIRECTORY)
  get_filename_component(_fjcontrib_prefix ${FastJetContrib_INCLUDE_DIR} DIRECTORY)
  foreach(_h ${_fjcontrib_headers})
    if(NOT EXISTS ${FastJetContrib_INCLUDE_DIR}/fastjet/contrib/${_h})
      message(WARNING "FastJet contrib in ${FastJetContrib_INCLUDE_DIR} has no ${_h}")
      set(FastJetContrib_INCLUDE_DIR FastJetContrib_INCLUDE_DIR-NOTFOUND)
    endif()
  endforeach()
endif()

find_library(FastJetContrib_fragile_LIBRARY fastjetcontribfragile
  HINTS ${_fjcontrib_prefix} ${_fjcontrib_hints} PATH_SUFFIXES lib lib64)
if(FastJetContrib_fragile_LIBRARY)
  set(_fjcontrib_lib_vars FastJetContrib_fragile_LIBRARY)
else()
  set(_fjcontrib_lib_vars)
  foreach(_lib Nsubjettiness ValenciaPlugin RecursiveTools)
    find_library(FastJetContrib_${_lib}_LIBRARY ${_lib}
      HINTS ${_fjcontrib_prefix} ${_fjcontrib_hints} PATH_SUFFIXES lib lib64)
    list(APPEND _fjcontrib_lib_vars FastJetContrib_${_lib}_LIBRARY)
  endforeach()
endif()

set(FastJetContrib_LIBRARIES)
foreach(_var ${_fjcontrib_lib_vars})
  if(${_var})
    get_filename_component(_real ${${_var}} REALPATH)
    list(APPEND FastJetContrib_LIBRARIES ${_real})
  endif()
endforeach()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FastJetContrib
  REQUIRED_VARS FastJetContrib_INCLUDE_DIR ${_fjcontrib_lib_vars})

if(FastJetContrib_FOUND)
  set(FastJetContrib_INCLUDE_DIRS ${FastJetContrib_INCLUDE_DIR})
  if(NOT TARGET FastJet::Contrib)
    add_library(FastJet::Contrib INTERFACE IMPORTED)
    set_target_properties(FastJet::Contrib PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${FastJetContrib_INCLUDE_DIRS}"
      INTERFACE_LINK_LIBRARIES "${FastJetContrib_LIBRARIES};FastJet::FastJet")
  endif()
endif()

mark_as_advanced(_fjcontrib_header_dir FastJetContrib_fragile_LIBRARY
  FastJetContrib_Nsubjettiness_LIBRARY FastJetContrib_ValenciaPlugin_LIBRARY FastJetContrib_RecursiveTools_LIBRARY)
