include_directories(
  ${CMAKE_SOURCE_DIR}
  ${ROOT_INCLUDE_DIRS}
  ${DelphesExternals_INCLUDE_DIR}
)

file(GLOB sources *.cc)
file(GLOB headers *.h)

add_library(modules OBJECT ${sources})

# install headers
install(FILES ${headers} DESTINATION include/modules)
