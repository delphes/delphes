find_package(pybind11)
if(NOT pybind11_FOUND)
  return()
endif()

file(GLOB sources *.cc)
file(GLOB headers *.h)

pybind11_add_module(DelphesPython ${sources})
target_link_libraries(DelphesPython PUBLIC Delphes)
target_include_directories(DelphesPython PRIVATE ${CMAKE_SOURCE_DIR} ${pybind11_INCLUDE_DIRS}
  ${ROOT_INCLUDE_DIRS} ${DelphesExternals_INCLUDE_DIR} # to be removed
)
set_target_properties(DelphesPython PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# install headers
install(FILES ${headers} DESTINATION include/DelphesPython COMPONENT python)
