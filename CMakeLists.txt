cmake_minimum_required(VERSION 3.20)
project(traminex_kernels LANGUAGES NONE)

# The default package is pure Python. The supported compiled profile is an
# explicit pybind11 build, so compiler discovery must stay inside this block.
option(TRAMINEX_BUILD_PYBIND11 "Build pybind11 bindings" OFF)

if (TRAMINEX_BUILD_PYBIND11)
  enable_language(CXX)
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)

  find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
  find_package(pybind11 CONFIG REQUIRED)

  # OpenMP parallelises the R-source-derived accel kernels. Keep compiler portability by
  # allowing a single-threaded build; the maintained Ubuntu profile verifies it.
  find_package(OpenMP)

  pybind11_add_module(_pb_kernels cpp/bindings/pybind11_module.cpp)
  target_include_directories(_pb_kernels PRIVATE cpp/src)
  target_compile_definitions(_pb_kernels PRIVATE PYBIND11_DETAILED_ERROR_MESSAGES)
  if (OpenMP_CXX_FOUND)
    target_link_libraries(_pb_kernels PRIVATE OpenMP::OpenMP_CXX)
    message(STATUS "_pb_kernels: OpenMP enabled (${OpenMP_CXX_FLAGS})")
  else()
    message(STATUS "_pb_kernels: OpenMP not found -- accel kernels build single-threaded")
  endif()

  install(TARGETS _pb_kernels LIBRARY DESTINATION traminex)
endif()
