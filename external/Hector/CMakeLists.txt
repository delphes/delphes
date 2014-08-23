include_directories(
  ${CMAKE_SOURCE_DIR}/external
  ${ROOT_INCLUDE_DIRS}
)

file(GLOB sources *.cc)

add_library(Hector OBJECT ${sources})
