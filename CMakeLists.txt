cmake_minimum_required(VERSION 3.14)
project(pmv_calc)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 REQUIRED)
find_package(CGAL REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(Boost REQUIRED)
find_package(OpenMP REQUIRED)

find_library(GMP_LIBRARY NAMES gmp REQUIRED)
find_library(MPFR_LIBRARY NAMES mpfr)

pybind11_add_module(pmv_calc 
    src/module.cpp 
    src/module_monte.cpp
)

include_directories(${Eigen3_INCLUDE_DIRS})

target_link_libraries(pmv_calc PRIVATE
    CGAL::CGAL
    Eigen3::Eigen
    Boost::boost
    OpenMP::OpenMP_CXX
    ${GMP_LIBRARY}
    ${MPFR_LIBRARY}
)
