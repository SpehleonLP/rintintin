cmake_minimum_required(VERSION 3.10)

# Set the project name and version
project(rintintin VERSION 1.0.0)


# Set C standard
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Same for C compiler
set(CMAKE_C_FLAGS_DEBUG "-g -O0 -DDEBUG")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g -DNDEBUG")
set(CMAKE_C_FLAGS_MINSIZEREL "-Os -DNDEBUG")

# Set output directory for the library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build)

# Collect all C source files from src/ directory
file(GLOB_RECURSE SRC_FILES "src/*.c")
file(GLOB_RECURSE SRC_HEADER_FILES "src/*.h")

# Collect all C source files from include/ directory (if any)
file(GLOB_RECURSE INCLUDE_C_FILES "include/*.h")

# Combine all source files
set(ALL_SOURCES ${SRC_FILES} ${SRC_HEADER_FILES} ${INCLUDE_C_FILES})

# Create the static library
add_library(rintintin STATIC ${ALL_SOURCES})

# Add include directory to the include path
target_include_directories(rintintin 
    PUBLIC 
        ${CMAKE_CURRENT_SOURCE_DIR}/include
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# Set the output name of the library
set_target_properties(rintintin PROPERTIES
    OUTPUT_NAME "rintintin"
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build/${CMAKE_BUILD_TYPE}
)

# Add debug postfix for debug builds
set_target_properties(rintintin PROPERTIES
    DEBUG_POSTFIX "_d"
)
