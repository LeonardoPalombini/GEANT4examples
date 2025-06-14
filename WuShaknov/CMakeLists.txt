# CMakeList.txt: progetto CMake per WS. Includere l'origine e definire
# qui la logica specifica del progetto.
#
cmake_minimum_required(VERSION 3.16...3.21)
project(WS)

option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/src)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(ws WS.cpp ${sources} ${headers})
target_link_libraries(ws ${Geant4_LIBRARIES})

#file(GLOB DATA_FILES "*.dat")
#file(COPY ${DATA_FILES} DESTINATION ${PROJECT_BINARY_DIR})

add_custom_target(WS DEPENDS ws)

install(TARGETS ws DESTINATION bin)


# Aggiungere l'origine all'eseguibile di questo progetto.
#add_executable (WS "WS.cpp" "WS.h")

#if (CMAKE_VERSION VERSION_GREATER 3.12)
 # set_property(TARGET WS PROPERTY CXX_STANDARD 20)
#endif()

# TODO: Aggiungere i test e, se necessario, installare le destinazioni.
