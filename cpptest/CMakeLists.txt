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


# IDK. maybe just glob all the cpps in vegafem
# add it to an uncompiled library and link that
# file(GLOB_RECURSE VegaFEM_SRCS "../VegaFEM-v4.0/libraries/*.cpp")
# # file(GLOB_RECURSE VegaFEM_HDRS "../VegaFEM-v4.0/libraries/*.h")
# message(${VegaFEM_SRCS})
# add_library(VegaFEMlib ${VegaFEM_SRCS})
# target_include_directories(VegaFEMlib PUBLIC  "../VegaFEM-v4.0/libraries/")
# target_link_libraries(lib PUBLIC VegaFEMlib)
# NOTE ALSO DOESNT WORK PROBABLY BECAUSE DAMN QUOTATOIN INCLUDES WITH OUT EVEN RELATIVE PATHS..

file(GLOB VegaFEM_libraries "${CMAKE_CURRENT_SOURCE_DIR}/../VegaFEM-v4.0/libraries/lib/*.a")
# message(${VegaFEM_libraries})
set(VegaFEM_include_directory "${CMAKE_CURRENT_SOURCE_DIR}/../VegaFEM-v4.0/libraries/include")
# target_link_libraries(lib PUBLIC ${VegaFEM_libraries})
# foreach(file ${VegaFEM_libraries})
#   message(${file})
#   target_link_libraries(lib PUBLIC ${file})
# endforeach()
# target_include_directories(lib PUBLIC ${VegaFEM_include_directory})
# include_directories(${VegaFEM_include_directory})

# APPROACH: COMBINE STATIC LIBRARIES INTO ONE
# add_custom_target(combined ALL COMMAND ${CMAKE_AR} rc libcombined.a $<TARGET_FILE:lib1> $<TARGET_FILE:lib2>)
# add_custom_target(VegaFEMc ALL COMMAND ${CMAKE_CXX_ARCHIVE_CREATE} rc ${CMAKE_STATIC_LIBRARY_PREFIX}VegaFEMcombined${CMAKE_STATIC_LIBRARY_SUFFIX} ${VegaFEM_libraries})
# add_custom_target(VegaFEMc ALL COMMAND ${CMAKE_CXX_ARCHIVE_CREATE} rc ${CMAKE_STATIC_LIBRARY_PREFIX}VegaFEMcombined${CMAKE_STATIC_LIBRARY_SUFFIX} ${VegaFEM_libraries})
# add_custom_command(
#     OUTPUT ${LIBNAME}
#     COMMAND lib.exe /OUT:${LIBNAME} $<TARGET_FILE:lib1> $<TARGET_FILE:lib2>
#     DEPENDS lib1 lib2
#     COMMENT "Combining libs..."
# )
# add_library(VegaFEM STATIC IMPORTED)
# set_property(TARGET VegaFEM PROPERTY IMPORTED_LOCATION VegaFEMc)
# target_link_libraries(lib PUBLIC VegaFEM)
# target_include_directories(lib PUBLIC ${VegaFEM_include_directory})


# stuff left to try: 
#  link individual libs into separate named targets
#  add dependencies to individual links based on errors

# set(VegaFEM_libnames mesh minivector volumetricMesh isotropicHyperelasticFEM)
# foreach(libname ${VegaFEM_libnames})
#   message(${libname})
#   file(GLOB ${libname}_SRCS "../VegaFEM-v4.0/libraries/${libname}/*.cpp")
#   # file(GLOB ${libname}_HDRS "../VegaFEM-v4.0/libraries/${libname}/*.h")
#   add_library(VegaFEM_${libname} ${${libname}_SRCS})
#   # set_target_properties(${VegaFEM_${libname}} PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
#   target_link_libraries(lib PUBLIC VegaFEM_${libname})
#   # target_include_directories(VegaFEM_${libname} PUBLIC "../VegaFEM-v4.0/libraries/${libname}")
#   # target_include_directories(lib PUBLIC "../VegaFEM-v4.0/libraries/${libname}")
#   target_include_directories(VegaFEM_${libname} PUBLIC "../VegaFEM-v4.0/libraries/${libname}")
# endforeach()



# # import, link and include VegaFEM
# set(VegaFEM_libnames mesh minivector volumetricMesh isotropicHyperelasticFEM)
# foreach(libname ${VegaFEM_libnames})
#   add_library(VegaFEM_${libname} STATIC IMPORTED)
#   set_target_properties(VegaFEM_${libname} PROPERTIES
#     IMPORTED_LOCATION "/home/gsperl/Projects/hFEMs/VegaFEM-v4.0/libraries/lib/lib${libname}.a"
#     INTERFACE_INCLUDE_DIRECTORIES "/home/gsperl/Projects/hFEMs/VegaFEM-v4.0/libraries/include"
#   )
# endforeach()
# set_target_properties(VegaFEM_isotropicHyperelasticFEM PROPERTIES INTERFACE_LINK_LIBRARIES VegaFEM_volumetricMesh)
# foreach(libname ${VegaFEM_libnames})
#   target_link_libraries(lib PUBLIC VegaFEM_${libname}) # or maybe only isotrop. bc interface
# endforeach()
# target_include_directories(lib PUBLIC "/home/gsperl/Projects/hFEMs/VegaFEM-v4.0/libraries/include")

# add_library(vega::isotropicHyperelasticFEM STATIC IMPORTED)
# set_target_properties(vega::isotropicHyperelasticFEM PROPERTIES
#   IMPORTED_LOCATION "/home/gsperl/Projects/hFEMs/VegaFEM-v4.0/libraries/lib/libisotropicHyperelasticFEM.a"
#   INTERFACE_INCLUDE_DIRECTORIES "/home/gsperl/Projects/hFEMs/VegaFEM-v4.0/libraries/include"
#   INTERFACE_LINK_LIBRARIES "/home/gsperl/Projects/hFEMs/VegaFEM-v4.0/libraries/lib/libvolumetricMesh.a"
# )
# target_link_libraries(lib PUBLIC vega::isotropicHyperelasticFEM)
# target_include_directories(lib PUBLIC "/home/gsperl/Projects/hFEMs/VegaFEM-v4.0/libraries/include")


# add_library(VegaFEM STATIC IMPORTED)
# set_target_properties(VegaFEM PROPERTIES
#   IMPORTED_LOCATION "${VegaFEM_libraries}"
#   INTERFACE_INCLUDE_DIRECTORIES "/home/gsperl/Projects/hFEMs/VegaFEM-v4.0/libraries/include"
#   # INTERFACE_INCLUDE_DIRECTORIES "${VegaFEM_include_directory}"
# )
# target_link_libraries(lib PUBLIC VegaFEM)



set(vega_dir ${CMAKE_CURRENT_SOURCE_DIR}/../vega4)
file(GLOB_RECURSE vega_SRCS "${vega_dir}/*.cpp")
add_library(vega ${vega_SRCS})
target_link_libraries(lib PUBLIC vega)
target_include_directories(lib PUBLIC ${vega_dir})

##------------ MAIN --------------##

add_executable(cpptest "main.cpp")
target_link_libraries(cpptest PRIVATE lib)

