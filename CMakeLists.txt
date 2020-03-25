CMAKE_MINIMUM_REQUIRED(VERSION 3.10)

PROJECT(OpenGL)

SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)
SET(THIRD_PARTY_PATH ${CMAKE_SOURCE_DIR}/libs)


################################### PROJECT ####################################

SET(TARGET_NAME ga_bounce)
SET(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include/)

FILE(GLOB_RECURSE SOURCE_FILES ${CMAKE_SOURCE_DIR}/src/*.cpp)
FILE(GLOB_RECURSE HEADER_FILES ${CMAKE_SOURCE_DIR}/include/*.hpp)
FILE(GLOB_RECURSE SHADER_FILES ${CMAKE_SOURCE_DIR}/shader/*.glsl)
FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/output)

INCLUDE_DIRECTORIES(
    ${INCLUDE_DIR}
)

ADD_EXECUTABLE(${TARGET_NAME} ${SOURCE_FILES} ${HEADER_FILES} ${SHADER_FILES})


############################## 3RD PARTIES LIBS ################################

SET(OpenGL_GL_PREFERENCE GLVND)

# effolkronium / glm / lodepng
SET(LIBS_INCLUDE_DIR ${THIRD_PARTY_PATH}/include)

FIND_PACKAGE(SDL2 REQUIRED)
FIND_PACKAGE(OpenGL REQUIRED)
FIND_PACKAGE(GLEW REQUIRED)
FIND_PACKAGE(Eigen3 REQUIRED)
ADD_SUBDIRECTORY(${THIRD_PARTY_PATH}/c3ga)
ADD_SUBDIRECTORY(${THIRD_PARTY_PATH}/libcpuid)
ADD_SUBDIRECTORY(${THIRD_PARTY_PATH}/imgui)

##########################

INCLUDE_DIRECTORIES(
    SYSTEM
    ${SDL2_INCLUDE_DIRS}
    ${OPENGL_INCLUDE_DIR}
    ${GLEW_INCLUDE_DIR}
    ${CPUID_INCLUDE_DIRS}
    ${IMGUI_INCLUDE_DIRS}
    ${LIBS_INCLUDE_DIR}
    ${EIGEN3_INCLUDE_DIR}
    ${C3GA_INCLUDE_DIRS}
)

TARGET_LINK_LIBRARIES(
    ${TARGET_NAME}
    ${SDL2_LIBRARIES}
    ${OPENGL_LIBRARIES}
    ${GLEW_LIBRARY}
    ${C3GA_LIBRARIES}
    ${CPUID_LIBRARIES}
    ${IMGUI_LIBRARIES}
)

################################# COMPILATION ##################################

TARGET_COMPILE_OPTIONS(
    ${TARGET_NAME} PRIVATE
    -g
    -O3
    -std=c++17
)
# Adding as much warning as possible on GNU gcc/g++ and Clang
IF (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    TARGET_COMPILE_OPTIONS(
        ${TARGET_NAME} PRIVATE
        -Wall
        -Wextra
    )
ELSEIF (CMAKE_CXX_COMPILER_ID MATCHES "[Cc]lang")
    TARGET_COMPILE_OPTIONS(
        ${TARGET_NAME} PRIVATE
        -Weverything
        -Werror
        -Wfatal-errors
    )
    IF (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0.0")
    ELSE ()
    ENDIF ()
ENDIF ()


# Disable some warnings
TARGET_COMPILE_OPTIONS(
    ${TARGET_NAME} PRIVATE
    -Wno-unused-parameter
    -Wno-sign-compare
    -Wno-shadow
)
IF (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    IF (CMAKE_CXX_COMPILER_VERSION VERSION_LESS_EQUAL "8.0.0")
    ELSE ()
    ENDIF ()
ELSEIF (CMAKE_CXX_COMPILER_ID MATCHES "[Cc]lang")
    TARGET_COMPILE_OPTIONS(
        ${TARGET_NAME} PRIVATE
        -Wno-padded
        -Wno-c++98-compat
        -Wno-weak-vtables
        -Wno-unknown-pragmas
        -Wno-exit-time-destructors
        -Wno-sign-conversion
    )
    IF (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0.0")
    ELSE ()
    ENDIF ()
ENDIF ()
