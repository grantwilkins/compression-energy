cmake_minimum_required(VERSION 3.19)
project(MGARDX CXX)

# Find MGARD package
find_package(mgard REQUIRED)

# Set up PAPI
set(PAPI_ROOT_DIR "/jet/home/gwilkins/papi/src" CACHE PATH "Root directory of PAPI installation")

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

# Find MPI package
find_package(MPI REQUIRED)

# Find OpenMP package
find_package(OpenMP REQUIRED)

# Add the executables
add_executable(mgardx mgardx.cpp)
add_executable(mgardx_omp mgardx_omp.cpp)

# Include directories
target_include_directories(mgardx PRIVATE ${PAPI_INCLUDE_DIR})
target_include_directories(mgardx_omp PRIVATE ${PAPI_INCLUDE_DIR} ${MPI_CXX_INCLUDE_DIRS})

# Link libraries
target_link_libraries(mgardx 
    PRIVATE 
    mgard::mgard
    ${PAPI_LIBRARY}
)

target_link_libraries(mgardx_omp
    PRIVATE 
    mgard::mgard
    ${PAPI_LIBRARY}
    MPI::MPI_CXX
    OpenMP::OpenMP_CXX
)

# Set C++ standard (optional, but recommended)
set_target_properties(mgardx mgardx_omp PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED ON
)

# Set MPI compiler for mgardx_omp
set_target_properties(mgardx_omp PROPERTIES
    COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}"
    LINK_FLAGS "${MPI_CXX_LINK_FLAGS}"
)

# Enable OpenMP for mgardx_omp
target_compile_options(mgardx_omp PRIVATE ${OpenMP_CXX_FLAGS})

# Use mpicxx as the compiler for mgardx_omp
set_target_properties(mgardx_omp PROPERTIES
    CXX_COMPILER ${MPI_CXX_COMPILER}
)
