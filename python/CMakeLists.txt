cmake_minimum_required(VERSION 3.5)
project(hFEMs LANGUAGES CXX VERSION 0.1)

SET(CMAKE_CXX_STANDARD 14)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pthread")


##----------- CPP SOURCE LIB -------------##

file(GLOB_RECURSE SRCS "../src/*.hpp" "../src/*.cpp") # WHEN HAVING CPP SOURCES MIGHT NOT NEED HPP TODO
add_library(lib ${SRCS})
# target_compile_definitions(simulation_lib PUBLIC NO_PARALLEL)

# important to be linked with python module
set_target_properties(lib PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
# for header only stuff TODO
set_target_properties(lib PROPERTIES LINKER_LANGUAGE CXX)


##------------ DEPENDENCIES --------------##

find_package (Eigen3 CONFIG REQUIRED)
# set (COMMON_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR})
# target_include_directories(lib PUBLIC ${COMMON_INCLUDE_DIRS})
target_link_libraries(lib PUBLIC Eigen3::Eigen)




set(vega_dir ${CMAKE_CURRENT_SOURCE_DIR}/../vega4)
file(GLOB_RECURSE vega_SRCS "${vega_dir}/*.cpp")
add_library(vega ${vega_SRCS})
set_target_properties(vega PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
target_link_libraries(lib PUBLIC vega)
target_include_directories(lib PUBLIC ${vega_dir})

##------------ PYBIND --------------##

find_package(pybind11 CONFIG REQUIRED)

# bind the cpp sources with pybind11
pybind11_add_module(_solver src/bindings.cpp)

# link the source to the python module
target_link_libraries(_solver PRIVATE lib)

# current guess: libary name here, module name in bindings.cpp, and module name in setup.py have to match
