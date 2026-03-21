find_package(Python COMPONENTS Interpreter Development)
if(NOT Python_FOUND)
  return()
endif()
message(NOTICE "Python v${Python_VERSION} found")

file(GLOB sources *.cc)
file(GLOB headers *.h)

add_library(DelphesPython SHARED ${sources})
target_link_libraries(DelphesPython PUBLIC ${Python_LIBRARIES})
target_include_directories(DelphesPython PRIVATE ${CMAKE_SOURCE_DIR} ${Python_INCLUDE_DIRS})

# install headers
install(FILES ${headers} DESTINATION include/DelphesPython COMPONENT python)
