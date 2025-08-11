# math/CMakeLists.txt

# Create a static library target for the math module
add_library(math_module STATIC
    Vector3.cpp
)

# Public include directory for the math module headers
target_include_directories(math_module
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}  # So users of math_module can do #include "Vector3.h"
)

# Compiler warnings and flags for math module
target_compile_options(math_module
    PRIVATE
        -Wall -Wextra -Wpedantic
)
