find_package(Pythia8)
if(NOT Pythia8_FOUND)
  return()
endif()

file(GLOB sources *.cc)
file(GLOB headers *.h)

add_library(DelphesPythia SHARED ${sources})
target_link_libraries(DelphesPythia PUBLIC Delphes ${PYTHIA8_LIBRARIES} ${CMAKE_DL_LIBS})
target_include_directories(DelphesPythia PRIVATE ${CMAKE_SOURCE_DIR}
  ${DelphesExternals_INCLUDE_DIR} # to be removed
)
set_target_properties(DelphesPythia PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# install headers
install(FILES ${headers} DESTINATION include/DelphesPythia COMPONENT pythia8)
# install library
install(TARGETS DelphesPythia DESTINATION lib COMPONENT pythia8)

# build all executables and put them into bin/
add_executable(DelphesPythia8 DelphesPythia8.cpp)
target_link_libraries(DelphesPythia8 PUBLIC DelphesPythia)
target_include_directories(DelphesPythia8 PRIVATE ${CMAKE_SOURCE_DIR} ${DelphesExternals_INCLUDE_DIR})
set_target_properties(DelphesPythia8 PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/readers)
install(TARGETS DelphesPythia8 DESTINATION bin)
