cmake_minimum_required(VERSION 3.9)

project(c3ga CXX)
set(CMAKE_BUILD_TYPE Release)

# include eigen
find_package(Eigen3 REQUIRED)
if(${EIGEN3_FOUND})
else()
    include_directories("/usr/include/eigen3") # manually specify the include location
endif()

option(BUILD_PYTHON "Build Python library" OFF)
if (BUILD_PYTHON)
    find_package(pybind11 REQUIRED)
endif()

# compilation flags
if (MSVC)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O2 -std=c++14")
else()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O2 -std=c++14 -Wno-return-local-addr")
endif()
set(CMAKE_CXX_FLAGS_DEBUG_INIT "-g")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftime-report ") # comilation profiling:  -Q -ftime-report
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-elide-constructors") # force the compilor to use our move 'copy constructor'

# files to compile
set(source_files src/c3ga/Mvec.cpp)
file(GLOB_RECURSE header_files src/c3ga/*.hpp)

# generate a library
if (MSVC)
	add_library(c3ga STATIC ${source_files} ${header_files})
else()
	add_library(c3ga SHARED ${source_files} ${header_files})
endif()

if (BUILD_PYTHON)
    pybind11_add_module(c3ga_py
        src/c3ga/PythonBindings.cpp
    )
    target_link_libraries(c3ga_py PRIVATE 
        c3ga)
endif()



# include directory path
include_directories(src)
include_directories(${EIGEN3_INCLUDE_DIR})
#
## install lib
#install(FILES ${header_files} ${source_files} DESTINATION /usr/local/include/c3ga)
#install(TARGETS c3ga
#        RUNTIME DESTINATION /usr/local/bin
#        LIBRARY DESTINATION /usr/local/lib
#        ARCHIVE DESTINATION /usr/local/lib)
#

SET(C3GA_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src/ PARENT_SCOPE)
SET(C3GA_LIBRARIES c3ga PARENT_SCOPE)
