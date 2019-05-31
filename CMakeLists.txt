cmake_minimum_required(VERSION 3.0.0)
project(BDUMP VERSION 0.1.0)

include(CTest)
enable_testing()

#add_executable(BDUMP main.cpp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()
include(${Geant4_USE_FILE})
option(WITH_ANALYSIS_USE "Build example with analysis objects" ON)
if(WITH_ANALYSIS_USE)
add_definitions(-DANALYSIS_USE)
else()
endif()

include_directories(${PROJECT_SOURCE_DIR}/include  
                     ${Geant4_INCLUDE_DIR} )
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(${PROJECT_NAME} main.cpp ${sources} ${headers})
target_link_libraries(${PROJECT_NAME} ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build basic_radioprotection. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(simulation_SCRIPTS
   vis.mac physics.mac
  )

foreach(_script ${simulation_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

install(TARGETS ${PROJECT_NAME} DESTINATION bin)
