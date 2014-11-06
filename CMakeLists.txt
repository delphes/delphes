cmake_minimum_required(VERSION 2.8)

project(Delphes)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -DDROP_CGAL")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")

# Declare ROOT dependency
find_package(ROOT COMPONENTS EG Eve Geom Gui GuiHtml GenVector Hist Physics Matrix Graf RIO Tree Gpad RGL MathCore)
include(${ROOT_USE_FILE})

if(NOT DEFINED CMAKE_INSTALL_LIBDIR)
  set(CMAKE_INSTALL_LIBDIR "lib")
endif()

function(DELPHES_GENERATE_DICTIONARY dictionary)
  if(${ROOT_VERSION} LESS 6.0)
    ROOT_GENERATE_DICTIONARY(${dictionary} ${ARGN})
  else()
    ROOT_GENERATE_DICTIONARY(${dictionary} MODULE ${dictionary} ${ARGN})
  endif()
endfunction()

# Declare position of all other externals needed
set(DelphesExternals_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/external)

add_subdirectory(classes)
add_subdirectory(converters)
add_subdirectory(display)
add_subdirectory(examples)
add_subdirectory(external)
add_subdirectory(modules)
add_subdirectory(readers)

add_library(Delphes SHARED
  $<TARGET_OBJECTS:classes> 
  $<TARGET_OBJECTS:display>
  $<TARGET_OBJECTS:modules> 
  $<TARGET_OBJECTS:ExRootAnalysis>
  $<TARGET_OBJECTS:fastjet>
  $<TARGET_OBJECTS:tcl>
  $<TARGET_OBJECTS:Hector>
)

target_link_Libraries(Delphes ${ROOT_LIBRARIES} ${ROOT_COMPONENT_LIBRARIES})

install(TARGETS Delphes DESTINATION lib)
