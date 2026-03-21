find_package(Python COMPONENTS Interpreter Development)
if(NOT Python_FOUND)
  return()
endif()
message(NOTICE "Python v${Python_VERSION} found")

find_package(pybind11)
if(NOT pybind11_FOUND)
  return()
endif()

file(GLOB sources *.cc)
file(GLOB headers *.h)

pybind11_add_module(DelphesPython ${sources})
target_link_libraries(DelphesPython PUBLIC Delphes ${Python_LIBRARIES} ${pybind11_LIBRARIES})
target_include_directories(DelphesPython PRIVATE ${CMAKE_SOURCE_DIR} ${Python_INCLUDE_DIRS} ${pybind11_INCLUDE_DIRS}
  ${ROOT_INCLUDE_DIRS} ${DelphesExternals_INCLUDE_DIR} # to be removed
)
set_target_properties(DelphesPython PROPERTIES CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# install headers
install(FILES ${headers} DESTINATION include/DelphesPython COMPONENT python)
