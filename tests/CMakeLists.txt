#top level CMakeList
cmake_minimum_required(VERSION 3.10)

#project name
project(darcyTest)

set(projectName "darcyTester.exe")
set(fileName "darcyTester.cpp")

# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" .
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/googletest-download"
)
execute_process(COMMAND "${CMAKE_COMMAND}" --build .
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/googletest-download"
)

# Prevent GoogleTest from overriding our compiler/linker options
# when building with Visual Studio
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This adds the following targets:
# gtest, gtest_main, gmock and gmock_main
add_subdirectory("${CMAKE_BINARY_DIR}/googletest-src"
                 "${CMAKE_BINARY_DIR}/googletest-build"
)

# The gtest/gmock targets carry header search path dependencies
# automatically when using CMake 2.8.11 or later. Otherwise we
# have to add them here ourselves.
if(CMAKE_VERSION VERSION_LESS 2.8.11)
    include_directories("${gtest_SOURCE_DIR}/include"
                        "${gmock_SOURCE_DIR}/include"
    )
endif()

#set compiler and flags
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
#compiler is only set, because jegors default compiler is set to intel not g++, because CATS
set(CMAKE_CXX_COMPILER /usr/bin/g++)
add_compile_options(-Wall -Wextra -pedantic -Ofast -march=native)

add_executable(${projectName} ${fileName})

#define directories where to lok for include files
target_include_directories(${projectName} PUBLIC
                            "${PROJECT_BINARY_DIR}"
                            "${PROJECT_SOURCE_DIR}/../libdarcy/include"
                            "${PROJECT_SOURCE_DIR}/../libdarcy_apps/include"
                            )
                            
#find eigen and its libray to our                            
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries (${projectName} Eigen3::Eigen gtest_main)
