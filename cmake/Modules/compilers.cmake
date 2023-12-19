

set(WARN_FLAGS "-Wall -Wundef -Wextra -Wshadow -Wno-register")

option(DO_COVERAGE "Enable compiler flags for measuring test coverage" OFF)

foreach(LANG CXX C)

  # general flags for both compilers/languages:
  set(CMAKE_${LANG}_FLAGS "${CMAKE_${LANG}_FLAGS} -finline -fPIC ${WARN_FLAGS}")

  # only for intel compilers:
  if(CMAKE_${LANG}_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_${LANG}_FLAGS "${CMAKE_${LANG}_FLAGS} -static-intel")
  endif()

  # only for GNU compilers:
  if(CMAKE_${LANG}_COMPILER_ID STREQUAL "GNU")
   set(CMAKE_${LANG}_FLAGS "${CMAKE_${LANG}_FLAGS} -ffast-math -Wno-stringop-truncation")
  endif()

  if(DO_COVERAGE)
    set(CMAKE_${LANG}_FLAGS "${CMAKE_${LANG}_FLAGS} -fprofile-arcs -ftest-coverage")
  endif()
endforeach()
 
