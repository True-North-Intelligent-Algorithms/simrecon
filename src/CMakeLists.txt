cmake_minimum_required(VERSION 3.10)
project(YourProjectName)

# Set the C standard to C11
set(CMAKE_C_STANDARD 11)

# Set the C++ standard to C++11
set(CMAKE_CXX_STANDARD 11)

find_package(OpenMP)
message(STATUS "OpenMP FOUND: " ${OPENMP_FOUND})
if(OPENMP_FOUND)
  # nvcc and openmp with CImg aren't working at the moment
  # https://github.com/dtschump/CImg/issues/286
  add_definitions(-Dcimg_use_openmp=0)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS
      "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


if(WIN32)
    set(PLATFORM win64)
else()
    set(PLATFORM linux64)
endif(WIN32)

set(PRIISM_LIB_PATH "${CMAKE_SOURCE_DIR}/IVE/${PLATFORM}/LIB")
set(PRIISM_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/IVE/${PLATFORM}/INCLUDE")
include_directories(${PRIISM_INCLUDE_PATH})
link_directories(${PRIISM_LIB_PATH})

message(STATUS "PRIISM_LIB_PATH: ${PRIISM_LIB_PATH}")

find_library(
    IMLIB
    NAMES imlib libimlib
    PATHS ${PRIISM_LIB_PATH} REQUIRED)

find_library(
    IVELIB
    NAMES ive libive
    PATHS ${PRIISM_LIB_PATH} REQUIRED)

include_directories( ${FFTW_INCLUDE_DIRS} )
include_directories( ${FFTW2_INCLUDE_DIRS})

# Create an executable for sirecon
add_executable(sirecon sirecon.c drift_n_phasestep.c small_math_funcs.c helpers.c)
add_executable(radialft radialft.cpp)

target_link_libraries(
  sirecon
  ${FFTW_FLOAT_LIB}
  ${IMLIB}
  ${IVELIB}
  ${LAPACK_LIBRARIES}
  ${BLAS_LIBRARIES}
  ${OpenMP_C_FLAGS}
 )

find_package(TIFF)

 target_link_libraries(
    radialft
    ${TIFF_LIBRARIES}
    ${FFTW_FLOAT_LIB}
    ${IMLIB}
    ${IVELIB}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES}
    ${OpenMP_C_FLAGS}
 )


# Create an executable for radialft
#add_executable(radialft mainfile.c)
#target_link_libraries(radialft PRIVATE sfftw srfftw ${IWLIBS} -ltiff)

# Create an executable for otf2d
#add_executable(otf2d radialft2d_Hesper.cpp tiffhandle.cpp)
#target_link_libraries(otf2d PRIVATE sfftw srfftw ${IWLIBS} -ltiff)

# Create an executable for wiener2d
#add_executable(wiener2d wiener2d.c)
#target_link_libraries(wiener2d PRIVATE sfftw srfftw ${IWLIBS} -ltiff)