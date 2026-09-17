cmake_minimum_required(VERSION 3.15)
project(piva-online-cpp CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Look for dependencies in the active conda environment first
if(DEFINED ENV{CONDA_PREFIX})
  list(PREPEND CMAKE_PREFIX_PATH $ENV{CONDA_PREFIX})
endif()

find_path(FFTW_INCLUDE_DIR fftw3.h REQUIRED)
find_library(FFTW_LIBRARY fftw3 REQUIRED)
find_path(SNDFILE_INCLUDE_DIR sndfile.h REQUIRED)
find_library(SNDFILE_LIBRARY sndfile REQUIRED)

add_library(piva_online INTERFACE)
target_include_directories(piva_online INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/include ${FFTW_INCLUDE_DIR})
target_link_libraries(piva_online INTERFACE ${FFTW_LIBRARY})

# Optional: parallelize the frequency bins (used for 4+ channels)
find_package(OpenMP COMPONENTS CXX)
if(OpenMP_CXX_FOUND)
  target_link_libraries(piva_online INTERFACE OpenMP::OpenMP_CXX)
  message(STATUS "OpenMP enabled for the online algorithms")
else()
  message(STATUS "OpenMP not found: the online algorithms run single-threaded")
endif()

add_executable(verify_parallel tests/verify_parallel.cpp)
target_link_libraries(verify_parallel PRIVATE piva_online)
target_compile_options(verify_parallel PRIVATE -O3 -Wall -Wextra)

add_executable(separate_online src/separate_online.cpp)
target_include_directories(separate_online PRIVATE ${SNDFILE_INCLUDE_DIR})
target_link_libraries(separate_online PRIVATE piva_online ${SNDFILE_LIBRARY})
target_compile_options(separate_online PRIVATE -O3 -Wall -Wextra)

# Python module, used to check the C++ code against the Python implementation
option(BUILD_PYTHON "Build the piva_online_cpp Python module" ON)
if(BUILD_PYTHON)
  find_package(Python COMPONENTS Interpreter Development REQUIRED)
  execute_process(
    COMMAND ${Python_EXECUTABLE} -c "import pybind11; print(pybind11.get_cmake_dir())"
    OUTPUT_VARIABLE pybind11_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
  find_package(pybind11 CONFIG REQUIRED)
  pybind11_add_module(piva_online_cpp python/bindings.cpp)
  target_link_libraries(piva_online_cpp PRIVATE piva_online)
  target_compile_options(piva_online_cpp PRIVATE -O3)
endif()
