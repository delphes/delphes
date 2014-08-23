include_directories(
  ${CMAKE_SOURCE_DIR}
  ${ROOT_INCLUDE_DIRS}
  ${DelphesExternals_INCLUDE_DIR}
)

file(GLOB sources *.cc)
file(GLOB headers *.h)
list(REMOVE_ITEM headers ${CMAKE_CURRENT_SOURCE_DIR}/ModulesLinkDef.h)
list(REMOVE_ITEM headers ${CMAKE_CURRENT_SOURCE_DIR}/Pythia8LinkDef.h)

DELPHES_GENERATE_DICTIONARY(ModulesDict ${headers} LINKDEF ModulesLinkDef.h)

# TODO: implement switch to enable Pythia8 if present
list(REMOVE_ITEM sources ${CMAKE_CURRENT_SOURCE_DIR}/PileUpMergerPythia8.cc)

add_library(modules OBJECT ${sources} ModulesDict.cxx)
