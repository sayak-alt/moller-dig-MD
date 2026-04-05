cmake_minimum_required(VERSION 3.10)
project(mollerdig)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --------------------------------------------------
# ROOT
# --------------------------------------------------
find_package(ROOT REQUIRED COMPONENTS Tree RIO)
include(${ROOT_USE_FILE})

# --------------------------------------------------
# Geant4 (REQUIRED now)
# --------------------------------------------------
find_package(Geant4 REQUIRED)

# --------------------------------------------------
# REMOLL path
# --------------------------------------------------
if(NOT DEFINED REMOLL_DIR)
    message(FATAL_ERROR "Please set REMOLL_DIR")
endif()

# --------------------------------------------------
# Include directories
# --------------------------------------------------
include_directories(
    ${ROOT_INCLUDE_DIRS}
    ${Geant4_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/include
    ${REMOLL_DIR}/include
)

# --------------------------------------------------
# Link directories
# --------------------------------------------------
link_directories(
    ${REMOLL_DIR}/lib
    ${REMOLL_DIR}/lib64
)

# --------------------------------------------------
# Sources
# --------------------------------------------------
set(SOURCES
    src/mollerdig.cxx
    src/remoll_tree.cxx
    src/remollDigAuxi.cxx
    src/mollerDigOutData.cxx
    src/GEMPlane.cxx
    src/GEMDet.cxx
    src/GEMSimDigitizer.cxx
)

add_executable(mollerdig ${SOURCES})

# --------------------------------------------------
# Link libraries
# --------------------------------------------------
target_link_libraries(mollerdig
    ${ROOT_LIBRARIES}
    ${Geant4_LIBRARIES}
    ${REMOLL_DIR}/lib64/libremoll${CMAKE_SHARED_LIBRARY_SUFFIX}
)

# --------------------------------------------------
# RPATH
# --------------------------------------------------
set_target_properties(mollerdig PROPERTIES
    BUILD_WITH_INSTALL_RPATH TRUE
#    INSTALL_RPATH "${REMOLL_DIR}/lib64"
    INSTALL_RPATH "${REMOLL_DIR}/lib64;${Geant4_DIR}/../../..;${ROOT_LIBRARY_DIR}"

)
