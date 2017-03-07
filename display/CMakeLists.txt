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

# install headers
install(FILES ${headers} DESTINATION include/display)

# install pcms if they are created
if (NOT ${ROOT_VERSION} VERSION_LESS "6.0.0")
  install(FILES
      ${PROJECT_BINARY_DIR}/display/libDisplayDict_rdict.pcm
      DESTINATION lib)
endif()
