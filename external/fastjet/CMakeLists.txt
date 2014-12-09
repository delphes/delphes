include_directories(${CMAKE_SOURCE_DIR}/external)

file(GLOB sources *.cc plugins/*/*.cc contribs/*/*.cc tools/*.cc)

add_library(fastjet OBJECT ${sources})
