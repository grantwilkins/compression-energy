cmake_minimum_required(VERSION 3.19)
project(MGARDX CXX)

# Find MGARD package
find_package(mgard REQUIRED)

# Set up PAPI
set(PAPI_ROOT_DIR "/home1/10191/gfw/papi/src" CACHE PATH "Root directory of PAPI installation")

find_path(PAPI_INCLUDE_DIR
    NAMES papi.h
    PATHS ${PAPI_ROOT_DIR}/include
    NO_DEFAULT_PATH
)

find_library(PAPI_LIBRARY
    NAMES papi
    PATHS ${PAPI_ROOT_DIR}/lib
    NO_DEFAULT_PATH
)

if(NOT PAPI_INCLUDE_DIR OR NOT PAPI_LIBRARY)
    message(FATAL_ERROR "PAPI not found. Please specify PAPI_ROOT_DIR.")
endif()

# Find HDF5
if(DEFINED ENV{TACC_HDF5_INC} AND DEFINED ENV{TACC_HDF5_LIB})
    set(HDF5_INCLUDE_DIR $ENV{TACC_HDF5_INC})
    set(HDF5_LIBRARY_DIR $ENV{TACC_HDF5_LIB})
    find_library(HDF5_LIBRARIES NAMES hdf5 PATHS ${HDF5_LIBRARY_DIR} NO_DEFAULT_PATH)
else()
    find_package(HDF5 REQUIRED)
endif()

# Find NetCDF
find_package(NetCDF REQUIRED)

# Add the executables
add_executable(mgard_x_io_test mgard_x_io_test.cpp)

# Include directories
target_include_directories(mgard_x_io_test PRIVATE 
    ${PAPI_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}
    ${NETCDF_INCLUDE_DIR}
)

# Link libraries
target_link_libraries(mgard_x_io_test 
    PRIVATE 
    mgard::mgard
    ${PAPI_LIBRARY}
    ${HDF5_LIBRARIES}
    ${NETCDF_LIBRARIES}
)

# Set C++ standard
set_target_properties(mgard_x_io_test PROPERTIES
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED ON
)

# Add OpenMP support
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(mgard_x_io_test PRIVATE OpenMP::OpenMP_CXX)
endif()

# Add MPI support
find_package(MPI REQUIRED)
target_link_libraries(mgard_x_io_test PRIVATE MPI::MPI_CXX)

# Add any additional flags
target_compile_options(mgard_x_io_test PRIVATE -D_POSIX_C_SOURCE=199309L)