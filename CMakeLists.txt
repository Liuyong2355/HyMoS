cmake_minimum_required(VERSION 3.15)
project(HyMoS VERSION 1.0.0 LANGUAGES CXX)

# ===== 编译选项 =====
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -march=x86-64 -fopenmp -DHAVE_INLINE")

# ===== 依赖库 =====
find_package(Boost REQUIRED COMPONENTS serialization system thread)
find_package(Eigen3 REQUIRED)
find_package(MPI QUIET)

if(NOT MPI_CXX_FOUND)
    message(STATUS "CMake FindMPI failed, falling back to mpicxx -show")
    find_program(MPI_CXX_COMPILER NAMES mpicxx HINTS /usr/bin /usr/lib/x86_64-linux-gnu/openmpi/bin)
    if(NOT MPI_CXX_COMPILER)
        message(FATAL_ERROR "MPI C++ compiler not found")
    endif()

    execute_process(
        COMMAND ${MPI_CXX_COMPILER} -show
        OUTPUT_VARIABLE MPI_CXX_SHOW
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    string(REGEX MATCHALL "(-I[^ ]+)|(-L[^ ]+)|(-l[^ ]+)" MPI_CXX_SHOW_FLAGS "${MPI_CXX_SHOW}")
    set(MPI_CXX_INCLUDE_DIRS "")
    set(MPI_CXX_LIBRARY_DIRS "")
    set(MPI_CXX_LIBRARIES "")

    foreach(flag IN LISTS MPI_CXX_SHOW_FLAGS)
        if(flag MATCHES "^-I(.+)")
            list(APPEND MPI_CXX_INCLUDE_DIRS "${CMAKE_MATCH_1}")
        elseif(flag MATCHES "^-L(.+)")
            list(APPEND MPI_CXX_LIBRARY_DIRS "${CMAKE_MATCH_1}")
        elseif(flag MATCHES "^-l(.+)")
            list(APPEND MPI_CXX_LIBRARIES "${flag}")
        endif()
    endforeach()

    find_path(MPI_CXX_HEADER_DIR mpi.h
        HINTS
            ${MPI_CXX_INCLUDE_DIRS}
            /usr/lib/x86_64-linux-gnu/openmpi/include
    )
    if(NOT MPI_CXX_HEADER_DIR)
        message(FATAL_ERROR "mpi.h not found for MPI fallback configuration")
    endif()

    list(APPEND MPI_CXX_INCLUDE_DIRS "${MPI_CXX_HEADER_DIR}")
    list(REMOVE_DUPLICATES MPI_CXX_INCLUDE_DIRS)

    if(MPI_CXX_LIBRARY_DIRS)
        link_directories(${MPI_CXX_LIBRARY_DIRS})
    endif()

    set(MPI_CXX_FOUND TRUE)
endif()

if(NOT MPI_CXX_FOUND)
    message(FATAL_ERROR "Could NOT find MPI C++ support")
endif()



# ===== 链接库 =====
set(COMMON_LIBS
    ${MPI_CXX_LIBRARIES}
    Boost::serialization
    Boost::system
    Boost::thread
    fftw3
    tbb
    m
    dl
)

# ===== 公共源文件 =====
set(RUNTIME_COMMON_SRC
    src/common/chameleon.cpp
    src/common/ConfigFile.cpp
    src/common/Controller.cpp
)

set(ENGINE_SRC
    src/engine/CaseRegistry.cpp
    src/engine/ConfigSnapshot.cpp
    src/engine/BuiltinLaunchers1D.cpp
    src/engine/BuiltinLaunchers2D.cpp
    src/engine/Engine.cpp
)

add_library(hymos_core SHARED
    ${RUNTIME_COMMON_SRC}
    ${ENGINE_SRC}
)
add_library(HyMoS::core ALIAS hymos_core)

target_include_directories(hymos_core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/engine>
    PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src/common
    ${CMAKE_SOURCE_DIR}/src/engine
    ${CMAKE_SOURCE_DIR}/src/1d
    ${CMAKE_SOURCE_DIR}/src/2d
    ${CMAKE_SOURCE_DIR}/cases/1d
    ${CMAKE_SOURCE_DIR}/cases/2d
    ${CMAKE_SOURCE_DIR}/apps
    ${EIGEN3_INCLUDE_DIR}
    ${MPI_CXX_INCLUDE_DIRS}
)

target_link_libraries(hymos_core PRIVATE ${COMMON_LIBS})
set_target_properties(hymos_core PROPERTIES
    OUTPUT_NAME hymos
    VERSION ${PROJECT_VERSION}
    # CouetteConfig carries an internal legacy-distribution output path. Keep
    # this plugin from binding to an already-loaded ABI-1/2/3/4 libhymos.
    SOVERSION 5
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
)

add_executable(hymos apps/hymos.cpp)
target_link_libraries(hymos PRIVATE HyMoS::core)
target_compile_definitions(hymos PRIVATE HYMOS_VERSION="${PROJECT_VERSION}")
set_target_properties(hymos PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

# Opt-in Baltamatica 2025 Linux plugin. Worker threads use only libhymos.
option(HYMOS_BUILD_BALTAMATICA_PLUGIN "Build the Baltamatica plugin" OFF)
if(HYMOS_BUILD_BALTAMATICA_PLUGIN)
    find_package(Baltamatica 4.1 REQUIRED HINTS /opt/Baltamatica)
    set(HYMOS_PLUGIN_BUNDLE_DIR "${CMAKE_BINARY_DIR}/plugin/HyMoS")
    add_baltamatica_plugin(hymos_baltamatica_plugin plugins/baltamatica/hymos_plugin.cpp)
    target_link_libraries(hymos_baltamatica_plugin PRIVATE HyMoS::core)
    target_include_directories(hymos_baltamatica_plugin PRIVATE ${CMAKE_SOURCE_DIR}/src/engine)
    set_target_properties(hymos_baltamatica_plugin PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY "${HYMOS_PLUGIN_BUNDLE_DIR}"
        BUILD_WITH_INSTALL_RPATH TRUE
        BUILD_RPATH "\$ORIGIN/lib;$<TARGET_FILE_DIR:Baltamatica::bex>"
        INSTALL_RPATH "\$ORIGIN/lib;$<TARGET_FILE_DIR:Baltamatica::bex>")
    add_custom_target(hymos_baltamatica_bundle ALL
        COMMAND ${CMAKE_COMMAND}
            "-DHYMOS_SOURCE_DIR=${CMAKE_SOURCE_DIR}"
            "-DHYMOS_BUNDLE_DIR=${HYMOS_PLUGIN_BUNDLE_DIR}"
            "-DHYMOS_CORE_FILE=$<TARGET_FILE:hymos_core>"
            "-DHYMOS_CORE_SONAME=$<TARGET_SONAME_FILE_NAME:hymos_core>"
            "-DHYMOS_VERSION=${PROJECT_VERSION}"
            -P "${CMAKE_SOURCE_DIR}/scripts/assemble_plugin.cmake"
        DEPENDS hymos_baltamatica_plugin hymos_core
        VERBATIM
        COMMENT "Assembling the HyMoS Baltamatica plugin")
endif()

# Opt-in integration checks; no extra dependency for ordinary solver builds.
option(HYMOS_BUILD_TESTS "Build and register library integration checks" OFF)
if(HYMOS_BUILD_TESTS)
    enable_testing()
    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    add_executable(hymos_config_file_test
        tests/library/config_file.cpp
        src/common/ConfigFile.cpp
        src/common/chameleon.cpp)
    target_include_directories(hymos_config_file_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/common)
    add_test(NAME config_file_contract COMMAND hymos_config_file_test)
    add_executable(hymos_couette_memory_test tests/library/couette_memory.cpp)
    target_link_libraries(hymos_couette_memory_test PRIVATE HyMoS::core)
    add_test(NAME couette_memory_and_cli
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/library/check_couette.py
            $<TARGET_FILE:hymos_couette_memory_test>
            $<TARGET_FILE:hymos>
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/nodealii/couette
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/nodealii/couette/input.txt)
    set_tests_properties(couette_memory_and_cli PROPERTIES TIMEOUT 240)
endif()
