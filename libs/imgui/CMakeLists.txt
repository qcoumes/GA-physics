CMAKE_MINIMUM_REQUIRED(VERSION 3.10)

PROJECT(ImGUI)

FILE(GLOB_RECURSE IMGUI_SOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
FILE(GLOB_RECURSE IMGUI_HEADER_FILES ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h)

ADD_LIBRARY(imgui SHARED ${IMGUI_SOURCE_FILES} ${IMGUI_HEADER_FILES})

SET(IMGUI_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include/)
INCLUDE_DIRECTORIES(
    ${IMGUI_INCLUDE_DIRS}
)


SET(IMGUI_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include/ PARENT_SCOPE)
SET(IMGUI_LIBRARIES imgui PARENT_SCOPE)
