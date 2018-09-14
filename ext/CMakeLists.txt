include(ExternalProject)

# Niels Lohmann's json library (single-header version).

add_library(ext-json INTERFACE)
target_include_directories(ext-json INTERFACE json/single_include)

# TCLAP command line parsing libary (header-only).

add_library(ext-tclap INTERFACE)
target_include_directories(ext-tclap INTERFACE tclap/include)

# Google benchmark for microbenchmarks:

check_git_submodule(gbench_avail google-benchmark QUIET)
if(gbench_avail)
    # Set up google benchmark as an external project.

    set(gbench_src_dir "${CMAKE_CURRENT_SOURCE_DIR}/google-benchmark")
    set(gbench_install_dir "${CMAKE_CURRENT_BINARY_DIR}/google-benchmark")
    set(gbench_cmake_args
        "-DCMAKE_BUILD_TYPE=release"
        "-DCMAKE_INSTALL_PREFIX=${gbench_install_dir}"
        "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
        "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")

    ExternalProject_Add(gbench
        # Add dummy DOWNLOAD_COMMAND to stop ExternalProject_Add terminating CMake if the
        # git submodule had not been udpated.
        DOWNLOAD_COMMAND "${CMAKE_COMMAND}" -E echo "Warning: ${gbench_src_dir} empty or missing."
        SOURCE_DIR "${gbench_src_dir}"
        CMAKE_ARGS "${gbench_cmake_args}"
        INSTALL_DIR "${gbench_install_dir}"
    )
    set_target_properties(gbench PROPERTIES EXCLUDE_FROM_ALL TRUE)

else()
    add_error_target(gbench
        "Building Google benchmark library"
        "The git submodule for google benchmark is not available")
endif()

add_library(ext-benchmark INTERFACE)
add_dependencies(ext-benchmark gbench)
target_include_directories(ext-benchmark INTERFACE "${gbench_install_dir}/include")
target_link_libraries(ext-benchmark INTERFACE "${gbench_install_dir}/lib/libbenchmark.a")

# Libraries, cmake configs/modules for contra.
# Contra requires Catch2, conduit, [...] to be able to be found via find_package().

set(ext-cmake-module-path ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_BINARY_DIR}/share/cmake")
set(ext-cmake-prefix-path ${CMAKE_PREFIX_PATH})

# Yes, really.
string(REPLACE ";" ":" _ext-cmake-module-path-with-colon "${ext-cmake-module-path}")
string(REPLACE ";" ":" _ext-cmake-prefix-path-with-colon "${ext-cmake-prefix-path}")
set(def-cmake-paths-with-colon
     "-DCMAKE_MODULE_PATH=${_ext-cmake-module-path-with-colon}"
     "-DCMAKE_PREFIX_PATH=${_ext-cmake-prefix-path-with-colon}")

# Catch library for contra

find_package(Catch2 QUIET)
if(NOT Catch2_FOUND)
    check_git_submodule(catch_avail catch)
    if(catch_avail)
        set(catch_src_dir "${CMAKE_CURRENT_SOURCE_DIR}/catch")
        set(catch_install_dir "${CMAKE_CURRENT_BINARY_DIR}")
        set(catch_cmake_args
            ${def-cmake-paths-with-colon}
            "-DCATCH_BUILD_TESTING=OFF"
            "-DCATCH_BUILD_EXAMPLES=OFF"
            "-DCATCH_INSTALL_DOCS=OFF"
            "-DCATCH_INSTALL_HELPERS=OFF"
            "-DCMAKE_BUILD_TYPE=release"
            "-DCMAKE_INSTALL_PREFIX=${catch_install_dir}"
            "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")

        ExternalProject_Add(catch
            SOURCE_DIR "${catch_src_dir}"
            INSTALL_DIR "${catch_install_dir}"
            STEP_TARGETS build
            LIST_SEPARATOR :
            CMAKE_ARGS "${catch_cmake_args}"
        )
    else()
        add_error_target(catch
            "Building Catch2 library"
            "The git submodule for conduit is not available")

    endif()
endif()

# Conduit library for contra and nesci
# Note: Findconduit.cmake provided in toplevel cmake subdirectory.

find_package(conduit QUIET)
if(NOT conduit_FOUND)
    check_git_submodule(conduit_avail conduit)
    if(conduit_avail)
        set(conduit_src_dir "${CMAKE_CURRENT_SOURCE_DIR}/conduit/src")
        set(conduit_install_dir "${CMAKE_CURRENT_BINARY_DIR}")
        set(conduit_cmake_args
            ${def-cmake-paths-with-colon}
            "-DENABLE_TESTS=OFF"
            "-DENABLE_DOCS=OFF"
            "-DCMAKE_SKIP_RPATH=ON"
            "-DCMAKE_BUILD_TYPE=release"
            "-DCMAKE_INSTALL_PREFIX=${conduit_install_dir}"
            "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
            "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")

        ExternalProject_Add(conduit
            SOURCE_DIR "${conduit_src_dir}"
            INSTALL_DIR "${conduit_install_dir}"
            STEP_TARGETS build
            LIST_SEPARATOR :
            CMAKE_ARGS "${conduit_cmake_args}"
        )
    else()
        add_error_target(conduit
            "Building conduit library"
            "The git submodule for conduit is not available")
    endif()
endif()

# ZeroMQ C++ binding for contra.

find_package(cppzmq QUIET)
if(NOT cppzmq_FOUND)
    check_git_submodule(cppzmq_avail cppzmq)
    if(cppzmq_avail)
        set(cppzmq_src_dir "${CMAKE_CURRENT_SOURCE_DIR}/cppzmq")
        set(cppzmq_install_dir "${CMAKE_CURRENT_BINARY_DIR}")
        set(cppzmq_cmake_args
            ${def-cmake-paths-with-colon}
            "-DENABLE_DRAFTS=OFF"
            "-DCPPZMQ_BUILD_TESTS=OFF"
            "-DCMAKE_BUILD_TYPE=release"
            "-DCMAKE_INSTALL_PREFIX=${cppzmq_install_dir}"
            "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
            "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")

        ExternalProject_Add(cppzmq
            SOURCE_DIR "${cppzmq_src_dir}"
            INSTALL_DIR "${cppzmq_install_dir}"
            STEP_TARGETS build
            LIST_SEPARATOR :
            CMAKE_ARGS "${cppzmq_cmake_args}"
        )
    else()
        add_error_target(cppzmq
            "Building cppzmq library"
            "The git submodule for cppzmq is not available")
    endif()
endif()

# Contra library

find_package(contra QUIET)
if(NOT contra_FOUND)
    check_git_submodule(contra_avail contra)
    if(contra_avail)
        set(contra_src_dir "${CMAKE_CURRENT_SOURCE_DIR}/contra")
        set(contra_install_dir "${CMAKE_CURRENT_BINARY_DIR}")
        set(contra_cmake_args
            ${def-cmake-paths-with-colon}
            "-DCONDUIT_DIR=${conduit_install_dir}"
            "-DCMAKE_BUILD_TYPE=release"
            "-DCMAKE_INSTALL_PREFIX=${contra_install_dir}"
            "-DCMAKE_PREFIX_PATH=${CMAKE_CURRENT_BINARY_DIR}"
            "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
            "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")

        ExternalProject_Add(contra
            SOURCE_DIR "${contra_src_dir}"
            INSTALL_DIR "${contra_install_dir}"
            DEPENDS catch conduit
            STEP_TARGETS build
            LIST_SEPARATOR :
            CMAKE_ARGS "${contra_cmake_args}"
        )

        # TODO: make ext-contra interface library here that depends upon contra target
    else()
        add_error_target(conta
            "Building contra library"
            "The git submodule for contra is not available")
    endif()
else()
    # TODO: How does contra export itself?
    add_library(ext-contra ALIAS contra::contra)
endif()

# Sphinx RTD theme

check_git_submodule(rtdtheme_avail sphinx_rtd_theme)
add_target_if(rtdtheme_avail
    check-sphinx_rtd_theme
    "Checking for Read the Docs Sphinx theme"
    "The git submodule for read the docs is not available")

add_library(ext-sphinx_rtd_theme INTERFACE)
add_dependencies(ext-sphinx_rtd_theme check-sphinx_rtd_theme)
target_include_directories(ext-sphinx_rtd_theme INTERFACE sphinx_rtd_theme)
