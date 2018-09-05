# -------------------------------------------------------
# util cmake build script for sgffabric.
# Written by: Salvatore
# Modified by: Salvatore
# -------------------------------------------------------

# -------------------------------------------------------
# Source directories containing all the necessary .cpp files
# -------------------------------------------------------
set(SYS_DIR ${CMAKE_SOURCE_DIR}/Source/sys)
set(SYS_SRC
${SYS_DIR}/SGF_System.cpp
${SYS_DIR}/linux/SGF_LinShared.cpp
${SYS_DIR}/linux/SGF_LinCpu.cpp
${SYS_DIR}/linux/SGF_Fpu.cpp
)

# -------------------------------------------------------
# Include directory
# -------------------------------------------------------

# -------------------------------------------------------
# module
# -------------------------------------------------------
#${SYS_SRC}add_library (system_module ${SYS_SRC})
