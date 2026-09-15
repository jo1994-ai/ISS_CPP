# Q15 (16-bit fixed-point) FFT example. No external libraries.
# Windows: run build.bat, or open this folder in Visual Studio.
cmake_minimum_required(VERSION 3.15)
project(fft_example_q15 CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_executable(fft_example_q15 fft_example_q15.cpp q15_fft.hpp wav16.hpp)
if(MSVC)
  target_compile_options(fft_example_q15 PRIVATE /W3 /EHsc)
else()
  target_compile_options(fft_example_q15 PRIVATE -Wall -Wextra)
endif()

# put the program in the project folder, next to audio/
set_target_properties(fft_example_q15 PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_SOURCE_DIR}
  RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_CURRENT_SOURCE_DIR})
set_property(TARGET fft_example_q15 PROPERTY VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT fft_example_q15)
