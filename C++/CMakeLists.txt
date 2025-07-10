cmake_minimum_required(VERSION 3.0.0)
project(pldrop VERSION 0.1.0)

set (CMAKE_CXX_STANDARD 11)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)


#set_target_properties(Eigen3::Eigen PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "C:/Users/pbarr/.vscode/Programs/Planetary-Supply-Drop/C++/lib/eigen-3.4.0")

include_directories("C:/Users/pbarr/.vscode/Programs/Planetary-Supply-Drop/C++/lib/eigen-3.3.7")
list(APPEND CMAKE_MODULE_PATH "C:/Users/pbarr/.vscode/Programs/Planetary-Supply-Drop/C++/lib/eigen-3.3.7/cmake")
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

get_target_property(EigenIncludeDirs Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
message(STATUS "Eigen include dirs: ${EigenIncludeDirs}")
set_target_properties(Eigen3::Eigen PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "C:/Users/pbarr/.vscode/Programs/Planetary-Supply-Drop/C++/lib/eigen-3.3.7")
 
include_directories("include/")
include_directories("lib/")

add_compile_options("/utf-8")

add_executable(pldrop
    src/main.cpp
    include/edl.h src/edl.cpp
    include/orbits.h src/orbits.cpp
    )

set(ENABLE_ECOS TRUE)
set(ENABLE_OSQP TRUE)
add_subdirectory(lib/Epigraph)

target_link_libraries(pldrop Eigen3::Eigen epigraph)

get_target_property(EigenIncludeDirs Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
message(STATUS "Eigen include dirs: ${EigenIncludeDirs}")
