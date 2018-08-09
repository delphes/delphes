include_directories(
  ${CMAKE_SOURCE_DIR}
  ${DelphesExternals_INCLUDE_DIR} 
)

file(GLOB executables RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

# TODO: implement switch to enable CMSSW, ProMC and ProIO if present
list(REMOVE_ITEM executables DelphesCMSFWLite.cpp)
list(REMOVE_ITEM executables DelphesProMC.cpp)
list(REMOVE_ITEM executables DelphesProIO.cpp)
if(NOT PYTHIA8_FOUND)
  list(REMOVE_ITEM executables DelphesPythia8.cpp)
endif()

# build all executables and put them into bin/
foreach(sourcefile ${executables})
  string(REPLACE ".cpp" "" name ${sourcefile})
  add_executable(${name} ${sourcefile})
  target_link_libraries(${name} Delphes)
  install(TARGETS ${name} DESTINATION bin)
endforeach()
