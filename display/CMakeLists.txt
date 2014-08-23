include_directories(
  ${CMAKE_SOURCE_DIR}
  ${ROOT_INCLUDE_DIRS}
  ${DelphesExternals_INCLUDE_DIR} 
)

file(GLOB sources *.cc)
file(GLOB headers *.h)
list(REMOVE_ITEM headers ${CMAKE_CURRENT_SOURCE_DIR}/DisplayLinkDef.h)

DELPHES_GENERATE_DICTIONARY(DisplayDict ${headers} LINKDEF DisplayLinkDef.h)

add_library(display OBJECT ${sources} DisplayDict.cxx)
