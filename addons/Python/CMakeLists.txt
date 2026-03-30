find_package(pybind11)
if(NOT pybind11_FOUND)
  return()
endif()

file(GLOB sources *.cc)
file(GLOB headers *.h)
file(GLOB examples examples/*.py)

pybind11_add_module(DelphesPython ${sources})
target_link_libraries(DelphesPython PRIVATE Delphes pybind11::module)

target_include_directories(DelphesPython PRIVATE ${CMAKE_SOURCE_DIR}
  ${DelphesExternals_INCLUDE_DIR} # to be removed
)
set_target_properties(DelphesPython PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# install headers
install(FILES ${headers} DESTINATION include/DelphesPython COMPONENT python)
# install library
install(TARGETS DelphesPython DESTINATION lib COMPONENT python)
# install examples
install(FILES ${examples} DESTINATION examples COMPONENT python)
