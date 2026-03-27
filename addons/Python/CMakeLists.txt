find_package(pybind11)
if(NOT pybind11_FOUND)
  return()
endif()

file(GLOB sources *.cc)
file(GLOB headers *.h)

pybind11_add_module(DelphesPython MODULE ${sources})
target_link_libraries(DelphesPython PUBLIC
  $<TARGET_OBJECTS:classes>
  $<TARGET_OBJECTS:modules>
  $<TARGET_OBJECTS:ExRootAnalysis>
  $<TARGET_OBJECTS:fastjet>
  $<TARGET_OBJECTS:tcl>
  $<TARGET_OBJECTS:Hector>
  $<TARGET_OBJECTS:PUPPI>
  $<TARGET_OBJECTS:TrackCovariance>
  ROOT::Core ROOT::Hist ROOT::Tree ROOT::MathCore ROOT::Physics ROOT::EG
  ROOT::ROOTNTuple
  ${YAML_CPP_LIBRARIES}
)
if(PYTHIA8_FOUND)
  target_link_libraries(DelphesPython PUBLIC ${PYTHIA8_LIBRARIES} ${CMAKE_DL_LIBS})
endif()

target_include_directories(DelphesPython PRIVATE ${CMAKE_SOURCE_DIR}
  ${DelphesExternals_INCLUDE_DIR} # to be removed
)
set_target_properties(DelphesPython PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# install headers
install(FILES ${headers} DESTINATION include/DelphesPython COMPONENT python)
# install library
install(TARGETS DelphesPython DESTINATION lib COMPONENT python)
