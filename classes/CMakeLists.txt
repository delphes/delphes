include_directories(
  ${CMAKE_SOURCE_DIR}
  ${ROOT_INCLUDE_DIRS}
  ${DelphesExternals_INCLUDE_DIR} 
)

file(GLOB sources *.cc)
file(GLOB headers *.h)
list(REMOVE_ITEM headers ${CMAKE_CURRENT_SOURCE_DIR}/ClassesLinkDef.h)

DELPHES_GENERATE_DICTIONARY(ClassesDict 
  ${CMAKE_CURRENT_SOURCE_DIR}/DelphesModule.h
  ${CMAKE_CURRENT_SOURCE_DIR}/DelphesFactory.h
  ${CMAKE_CURRENT_SOURCE_DIR}/SortableObject.h
  ${CMAKE_CURRENT_SOURCE_DIR}/DelphesClasses.h
  LINKDEF ClassesLinkDef.h
)

add_library(classes OBJECT ${sources} ClassesDict.cxx)

# install public headers
install(FILES ${headers} DESTINATION include/classes)

# install pcms if they are created
if (NOT ${ROOT_VERSION} VERSION_LESS "6.0.0")
  install(FILES
      ${PROJECT_BINARY_DIR}/classes/libClassesDict_rdict.pcm
      DESTINATION lib)
endif()
