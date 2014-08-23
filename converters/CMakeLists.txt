include_directories(
  ${CMAKE_SOURCE_DIR}
  ${DelphesExternals_INCLUDE_DIR} 
)

# add all executables as targets
file(GLOB executables RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

foreach(sourcefile ${executables})
  string(REPLACE ".cpp" "" name ${sourcefile})
  add_executable(${name} ${sourcefile})
  target_link_libraries(${name} Delphes)
  install(TARGETS ${name} DESTINATION bin)
endforeach()
