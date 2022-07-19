include_directories(${CMAKE_SOURCE_DIR}/external)

file(GLOB sources *.cc)
file(GLOB headers *.h)

# install headers
install(FILES ${headers} DESTINATION include/TrackCovariance)

add_library(TrackCovariance OBJECT ${sources})
