include_directories(${CMAKE_SOURCE_DIR}/external)

file(GLOB sources *.c)

add_library(tcl OBJECT ${sources})
