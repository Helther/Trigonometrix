cmake_minimum_required(VERSION 3.5)

set(project_name trigTest)
project(${project_name} LANGUAGES CXX)

add_library(project_options INTERFACE)
target_compile_features(project_options INTERFACE cxx_std_20)

#add_compile_options(-ffast-math -mfma) #TODO for later

add_executable(${project_name}
    src/main.cpp
    src/trigonometry.hpp
    src/lut_generator.hpp
    src/float_table.hpp
    src/double_table.hpp
    src/polynomials_coeffs.hpp
    )


target_link_libraries(${project_name} PRIVATE project_options)

option(ENABLE_TESTING "Enable Test Builds" ON)

if(ENABLE_TESTING)
  message(
    "Building Tests."
  )
  add_subdirectory(tests)
  add_test(SinCosTests tests/tests 1)
endif()
