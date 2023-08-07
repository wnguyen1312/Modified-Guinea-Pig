
file(GLOB MODULES  EDMwriterCPP.cc besselCPP.cc beamParametersCPP.cc beamCPP.cc backgroundCPP.cc 
pairsCPP.cc minijetsCPP.cc mathematicalTools.cc lumiCPP.cc 
jetParameterCPP.cc guineapigCPP.cc grvCPP.cc gridCPP.cc
fourierCPP.cc fieldCPP.cc bhabhaSamplesCPP.cc timerCPP.cc 
switchesCPP.cc splineCPP.cc rndmCPP.cc resultsCPP.cc 
principal.cc physicalTools.cc particlesCPP.cc particleBeamCPP.cc 
readData.cc memory.cc tridentCPP.cc option_args.cc)

# Set up header macro generations
include(guinea_config)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/config.h)

add_executable(guinea ${MODULES})

option(GUINEA_USE_EDM4HEP "Generate edm4hep output" OFF)

if(GUINEA_USE_EDM4HEP)
  find_package(EDM4HEP REQUIRED)   
  find_package(podio)
  target_link_libraries(guinea EDM4HEP::edm4hep EDM4HEP::utils podio::podioRootIO)
  target_compile_definitions(guinea PUBLIC USE_EDM4HEP)
endif()
	

if(FFTW2)
   if(FFTW2_FOUND)
      target_link_libraries(guinea ${FFTW2_LIBRARIES})
      include_directories(${FFTW2_INCLUDES})
   else()
      target_link_libraries(guinea fftw m)
   endif()
elseif(FFTW3)
   if(FFTW_FOUND)
      target_link_libraries(guinea ${FFTW_LIBRARIES})
      include_directories(${FFTW_INCLUDES})
   else()
      target_link_libraries(guinea fftw3 m)
   endif()
endif()

# Installing targets and headers:
install(TARGETS guinea
  ARCHIVE
  DESTINATION lib
  COMPONENT Libraries
  RUNTIME
  DESTINATION bin
  COMPONENT Runtime)

file(GLOB headers *.h)
install(FILES ${headers}
  DESTINATION include/guinea-pig
  COMPONENT Headers)
