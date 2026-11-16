# Find an installed FastJet (with the SISCone and CDFCones plugins)
#
# Hints: FASTJET_ROOT (CMake or environment variable), fastjet-config in PATH
#
# Defines
#   FastJet_FOUND, FastJet_VERSION, FastJet_INCLUDE_DIRS, FastJet_LIBRARIES
#   the imported target FastJet::FastJet
#
# The include directory and libraries are resolved to their real location, so
# that the include directory of a software stack view (which typically also
# contains the headers of an installed Delphes) is not put on the include path.

find_program(FASTJET_CONFIG_EXECUTABLE fastjet-config
  HINTS ${FASTJET_ROOT}/bin $ENV{FASTJET_ROOT}/bin)

set(_fastjet_hints ${FASTJET_ROOT} $ENV{FASTJET_ROOT})
if(FASTJET_CONFIG_EXECUTABLE)
  execute_process(COMMAND ${FASTJET_CONFIG_EXECUTABLE} --prefix
    OUTPUT_VARIABLE _fastjet_prefix OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  list(APPEND _fastjet_hints ${_fastjet_prefix})
endif()

find_path(_fastjet_header_dir fastjet/ClusterSequence.hh
  HINTS ${_fastjet_hints} PATH_SUFFIXES include)

if(_fastjet_header_dir)
  get_filename_component(_fastjet_header ${_fastjet_header_dir}/fastjet/ClusterSequence.hh REALPATH)
  get_filename_component(_fastjet_header_dir ${_fastjet_header} DIRECTORY)
  get_filename_component(FastJet_INCLUDE_DIR ${_fastjet_header_dir} DIRECTORY)
  get_filename_component(_fastjet_prefix ${FastJet_INCLUDE_DIR} DIRECTORY)

  set(_fastjet_config_h ${FastJet_INCLUDE_DIR}/fastjet/config_auto.h)
  if(EXISTS ${_fastjet_config_h})
    file(STRINGS ${_fastjet_config_h} _fastjet_version_line REGEX "^#define FASTJET_PACKAGE_VERSION ")
    string(REGEX REPLACE ".*\"(.*)\".*" "\\1" FastJet_VERSION "${_fastjet_version_line}")
    foreach(_plugin SISCONE CDFCONES)
      file(STRINGS ${_fastjet_config_h} _fastjet_plugin_line REGEX "^#define FASTJET_ENABLE_PLUGIN_${_plugin} ")
      if(NOT _fastjet_plugin_line)
        set(_fastjet_missing_plugins "${_fastjet_missing_plugins} ${_plugin}")
      endif()
    endforeach()
  endif()
endif()

set(FastJet_LIBRARIES)
set(_fastjet_lib_vars)
foreach(_lib fastjettools fastjetplugins siscone_spherical siscone fastjet)
  find_library(FastJet_${_lib}_LIBRARY ${_lib}
    HINTS ${_fastjet_prefix} ${_fastjet_hints} PATH_SUFFIXES lib lib64)
  if(FastJet_${_lib}_LIBRARY)
    get_filename_component(_real ${FastJet_${_lib}_LIBRARY} REALPATH)
    list(APPEND FastJet_LIBRARIES ${_real})
  endif()
  list(APPEND _fastjet_lib_vars FastJet_${_lib}_LIBRARY)
endforeach()

# Additional libraries the FastJet libraries depend on (e.g. libgfortran for the
# PxCone plugin), which are not necessarily recorded in the shared libraries
if(FASTJET_CONFIG_EXECUTABLE)
  execute_process(COMMAND ${FASTJET_CONFIG_EXECUTABLE} --libs --plugins
    OUTPUT_VARIABLE _fastjet_config_libs OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  separate_arguments(_fastjet_config_libs UNIX_COMMAND "${_fastjet_config_libs}")
  foreach(_arg ${_fastjet_config_libs})
    if(_arg MATCHES "^-l(.+)$")
      set(_name ${CMAKE_MATCH_1})
      if(NOT _name MATCHES "^(fastjet|fastjettools|fastjetplugins|siscone|siscone_spherical)$")
        list(APPEND FastJet_LIBRARIES ${_name})
      endif()
    endif()
  endforeach()
  list(REMOVE_DUPLICATES FastJet_LIBRARIES)
endif()

if(_fastjet_missing_plugins)
  message(WARNING "FastJet in ${FastJet_INCLUDE_DIR} was built without the plugins:${_fastjet_missing_plugins}")
  set(FastJet_INCLUDE_DIR FastJet_INCLUDE_DIR-NOTFOUND)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FastJet
  REQUIRED_VARS FastJet_INCLUDE_DIR ${_fastjet_lib_vars}
  VERSION_VAR FastJet_VERSION)

if(FastJet_FOUND)
  set(FastJet_INCLUDE_DIRS ${FastJet_INCLUDE_DIR})
  if(NOT TARGET FastJet::FastJet)
    add_library(FastJet::FastJet INTERFACE IMPORTED)
    set_target_properties(FastJet::FastJet PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${FastJet_INCLUDE_DIRS}"
      INTERFACE_LINK_LIBRARIES "${FastJet_LIBRARIES}")
  endif()
endif()

mark_as_advanced(FASTJET_CONFIG_EXECUTABLE _fastjet_header_dir ${_fastjet_lib_vars})
