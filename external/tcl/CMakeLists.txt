include_directories(${DelphesExternals_INCLUDE_DIR})

file(GLOB sources *.c)

add_library(tcl OBJECT ${sources})
