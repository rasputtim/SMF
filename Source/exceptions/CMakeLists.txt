# -------------------------------------------------------
# util cmake build script for sgffabric.
# Written by: Salvatore
# Modified by: Salvatore
# -------------------------------------------------------

# -------------------------------------------------------
# Source directories containing all the necessary .cpp files
# -------------------------------------------------------
set(EXCEPTIONS_DIR ${CMAKE_SOURCE_DIR}/Source/exceptions)
set(EXCEPTIONS_SRC
${EXCEPTIONS_DIR}/SGF_Error.cpp
${EXCEPTIONS_DIR}/SGF_Exception.cpp
${EXCEPTIONS_DIR}/SGF_LoadException.cpp
${EXCEPTIONS_DIR}/SGF_MugenException.cpp
${EXCEPTIONS_DIR}/SGF_NetException.cpp
${EXCEPTIONS_DIR}/SGF_TokenException.cpp
${EXCEPTIONS_DIR}/SGF_FileExceptions.cpp)
${EXCEPTIONS_DIR}/SGF_ScriptException.cpp)
# -------------------------------------------------------
# Include directory
# -------------------------------------------------------

# -------------------------------------------------------
# module
# -------------------------------------------------------
#add_library (exceptions_module ${EXCEPTIONS_SRC})
