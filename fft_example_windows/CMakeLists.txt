# Builds the bundled FFTW (third_party/fftw-3.3.11.tar.gz) as a static library
# and links it into fft_example. Nothing is downloaded.
#
# Visual Studio 2019/2022:  run build.bat, or File > Open > Folder...
cmake_minimum_required(VERSION 3.18)
project(fft_example C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

#### FFTW from the bundled source
set(FFTW_VERSION 3.3.11)
set(FFTW_SRC ${CMAKE_BINARY_DIR}/fftw-${FFTW_VERSION})
if(NOT EXISTS ${FFTW_SRC}/CMakeLists.txt)
  message(STATUS "Extracting third_party/fftw-${FFTW_VERSION}.tar.gz")
  file(ARCHIVE_EXTRACT
       INPUT ${CMAKE_CURRENT_SOURCE_DIR}/third_party/fftw-${FFTW_VERSION}.tar.gz
       DESTINATION ${CMAKE_BINARY_DIR})
endif()

# static double-precision library only
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(DISABLE_FORTRAN ON CACHE BOOL "" FORCE)

# SIMD on 64-bit Intel/AMD; FFTW checks at run time which ones the CPU has
string(TOLOWER "${CMAKE_SYSTEM_PROCESSOR}" ARCH)
if(ARCH MATCHES "amd64|x86_64" AND NOT CMAKE_GENERATOR_PLATFORM MATCHES "ARM")
  set(ENABLE_SSE2 ON CACHE BOOL "" FORCE)
  set(ENABLE_AVX ON CACHE BOOL "" FORCE)
  set(ENABLE_AVX2 ON CACHE BOOL "" FORCE)
endif()

# FFTW's CMakeLists asks for CMake 3.0, which CMake 4 no longer accepts
set(CMAKE_POLICY_VERSION_MINIMUM 3.5)

# FFTW's own /bigobj setting comes after its library target, so set it here
if(MSVC)
  add_compile_options($<$<COMPILE_LANGUAGE:C>:/bigobj>)
endif()

add_subdirectory(${FFTW_SRC} ${CMAKE_BINARY_DIR}/fftw-build EXCLUDE_FROM_ALL)

#### the example
add_executable(fft_example fft_example.cpp wav_io.hpp)
target_include_directories(fft_example PRIVATE ${FFTW_SRC}/api)
target_link_libraries(fft_example PRIVATE fftw3)
if(MSVC)
  target_compile_options(fft_example PRIVATE /W3 /EHsc /utf-8)
  target_compile_definitions(fft_example PRIVATE _CRT_SECURE_NO_WARNINGS)
else()
  target_compile_options(fft_example PRIVATE -Wall -Wextra)
endif()

# put fft_example.exe in the project folder, next to audio/
set_target_properties(fft_example PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_SOURCE_DIR}
  RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_CURRENT_SOURCE_DIR})

# run from the project folder inside Visual Studio as well
set_property(TARGET fft_example PROPERTY VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT fft_example)
