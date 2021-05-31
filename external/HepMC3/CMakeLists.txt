include_directories(
  ${CMAKE_SOURCE_DIR}/external
  ${ROOT_INCLUDE_DIRS}
)

file(GLOB src_sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)
if(${ROOT_VERSION} LESS 6.0)
  message(STATUS "ROOT was not found or the version of ROOT ->${ROOT_VERSION}<- is too old. HepMC3 will be built without root I/O")
  add_library(HepMC3 OBJECT ${src_sources})
  target_include_directories(HepMC3 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include) 
else()
  message(STATUS "ROOT was found. The version of ROOT is ->${ROOT_VERSION}<-. HepMC3 will be built with root I/O")
  file(GLOB rootIO_sources ${CMAKE_CURRENT_SOURCE_DIR}/rootIO/src/*.cc)
  ROOT_GENERATE_DICTIONARY(G__HepMC3rootIODelphes ${CMAKE_CURRENT_SOURCE_DIR}/rootIO/include/rootIO_Classes.hh LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/rootIO/include/LinkDef.hh NOINSTALL OPTIONS -inlineInputHeader)
  add_library(HepMC3 OBJECT ${src_sources} ${rootIO_sources} G__HepMC3rootIODelphes.cxx)
  target_include_directories(HepMC3 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/rootIO/include ${ROOT_INCLUDE_DIRS})
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libHepMC3rootIODelphes_rdict.pcm DESTINATION  ${CMAKE_INSTALL_LIBDIR})
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libHepMC3rootIODelphes.rootmap   DESTINATION  ${CMAKE_INSTALL_LIBDIR})
endif()
