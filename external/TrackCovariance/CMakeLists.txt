include_directories(
  ${DelphesExternals_INCLUDE_DIR}
  ${ROOT_INCLUDE_DIRS}
)

file(GLOB sources *.cc)
file(GLOB headers *.h)

# install headers
install(FILES ${headers} DESTINATION include/TrackCovariance)

add_library(TrackCovariance OBJECT ${sources})
