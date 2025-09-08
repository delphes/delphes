include_directories(BEFORE ${DelphesExternals_INCLUDE_DIR})

file(GLOB sources *.cc plugins/*/*.cc contribs/*/*.cc tools/*.cc)

add_library(fastjet OBJECT ${sources})
