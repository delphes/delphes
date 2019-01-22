include_directories(
  ${CMAKE_SOURCE_DIR}
  ${ROOT_INCLUDE_DIRS}
  ${DelphesExternals_INCLUDE_DIR}
)

file(GLOB sources *.cc)
file(GLOB headers *.h)
list(REMOVE_ITEM headers ${CMAKE_CURRENT_SOURCE_DIR}/FastJetLinkDef.h)
list(REMOVE_ITEM headers ${CMAKE_CURRENT_SOURCE_DIR}/ModulesLinkDef.h)
list(REMOVE_ITEM headers ${CMAKE_CURRENT_SOURCE_DIR}/Pythia8LinkDef.h)

DELPHES_GENERATE_DICTIONARY(FastJetDict ${headers} LINKDEF FastJetLinkDef.h)
DELPHES_GENERATE_DICTIONARY(ModulesDict ${headers} LINKDEF ModulesLinkDef.h)

if(PYTHIA8_FOUND)
  DELPHES_GENERATE_DICTIONARY(Pythia8Dict ${headers} LINKDEF Pythia8LinkDef.h)
  add_library(modules OBJECT ${sources} FastJetDict.cxx ModulesDict.cxx Pythia8Dict.cxx)
else()
  list(REMOVE_ITEM sources ${CMAKE_CURRENT_SOURCE_DIR}/PileUpMergerPythia8.cc)
  add_library(modules OBJECT ${sources} FastJetDict.cxx ModulesDict.cxx)
endif()

# install headers
install(FILES ${headers} DESTINATION include/modules)

# install pcms if they are created
if (NOT ${ROOT_VERSION} VERSION_LESS "6.0.0")
  install(FILES
    ${PROJECT_BINARY_DIR}/modules/libModulesDict_rdict.pcm
    ${PROJECT_BINARY_DIR}/modules/libFastJetDict_rdict.pcm
    DESTINATION lib)
  if(PYTHIA8_FOUND)
    install(FILES
      ${PROJECT_BINARY_DIR}/modules/libPythia8Dict_rdict.pcm
      DESTINATION lib)
  endif()
endif()
