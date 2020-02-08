SET(GLM_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/libs/include/)
IF (EXISTS ${GLM_INCLUDE_DIRS})

ELSE ()
    MESSAGE(FATAL_ERROR "include for lib glm invalid or not found: ${GLM_INCLUDE_DIRS}")
ENDIF ()
