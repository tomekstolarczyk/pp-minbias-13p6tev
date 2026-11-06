cmake_minimum_required(VERSION 3.18)
project(pp-minbias-13p6tev LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# find Pythia 8
find_program(PYTHIA8_CONFIG pythia8-config REQUIRED)
execute_process(COMMAND ${PYTHIA8_CONFIG} --includedir
    OUTPUT_VARIABLE PYTHIA8_INCLUDE_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${PYTHIA8_CONFIG} --libdir
    OUTPUT_VARIABLE PYTHIA8_LIB_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
find_library(PYTHIA8_LIBRARY pythia8 HINTS ${PYTHIA8_LIB_DIR} REQUIRED)

# find ROOT: RIO provides TFile, Tree provides TTree, TreePlayer provides TTreeReader, Hist provides TH1D
find_package(ROOT REQUIRED COMPONENTS RIO Tree TreePlayer Hist)

add_executable(generate src/generate.cc)
target_compile_options(generate PRIVATE -Wall -Wextra -Wpedantic)
target_include_directories(generate PRIVATE ${PYTHIA8_INCLUDE_DIR})
target_link_libraries(generate PRIVATE ${PYTHIA8_LIBRARY} ROOT::RIO ROOT::Tree)

add_executable(analysis src/analysis.cc)
target_compile_options(analysis PRIVATE -Wall -Wextra -Wpedantic)
target_link_libraries(analysis PRIVATE ROOT::RIO ROOT::Tree ROOT::TreePlayer ROOT::Hist)
