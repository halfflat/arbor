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

# Sphinx RTD theme

check_git_submodule(rtdtheme_avail sphinx_rtd_theme)
add_target_if(rtdtheme_avail
    check-sphinx_rtd_theme
    "Checking for Read the Docs Sphinx theme"
    "The git submodule for read the docs is not available")

add_library(ext-sphinx_rtd_theme INTERFACE)
add_dependencies(ext-sphinx_rtd_theme check-sphinx_rtd_theme)
target_include_directories(ext-sphinx_rtd_theme INTERFACE sphinx_rtd_theme)

# Libraries, cmake configs/modules for contra.
# Contra requires Catch2, conduit, [...] to be able to be found via find_package().

# Force re-configure of cmake if external packages are built.
set(reconfigure_tag "${CMAKE_CURRENT_BINARY_DIR}/_force_reconfigure.cmake")
if(NOT EXISTS "${reconfigure_tag}")
    file(TOUCH "${reconfigure_tag}")
endif()
include("${reconfigure_tag}")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_BINARY_DIR}/share/cmake")
list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_BINARY_DIR}")

include(CMakeParseArguments)
function(add_external_cmake_submodule sub)
    cmake_parse_arguments(arg "" "SOURCE_DIR;TARGET;UNLESS" "CMAKE_ARGS;DEPENDS" ${ARGN})

    set(target "${sub}")
    if(arg_TARGET)
        set(target "${arg_TARGET}")
    endif()

    if(${arg_UNLESS})
        # Make dummy target for building
        add_custom_target("${target}" COMMAND true COMMENT "")
        message(STATUS "${arg_UNLESS} true; making dummy target ${target}")
    else()
        message(STATUS "${arg_UNLESS} false; making external project target ${target}")
        set(_submodule_avail)
        check_git_submodule(_submodule_avail ${sub})
        if(_submodule_avail)
            set(cmake_module_path ${CMAKE_MODULE_PATH})
            set(cmake_prefix_path ${CMAKE_PREFIX_PATH})
            string(REPLACE ";" ":" cmake_module_path "${cmake_module_path}")
            string(REPLACE ";" ":" cmake_prefix_path "${cmake_prefix_path}")

            set(src_dir "${CMAKE_CURRENT_SOURCE_DIR}/${sub}")
            if(arg_SOURCE_DIR)
                set(src_dir "${arg_SOURCE_DIR}")
            endif()

            set(install_dir "${CMAKE_CURRENT_BINARY_DIR}")

            set(cmake_args
                "-DCMAKE_MODULE_PATH=${cmake_module_path}"
                "-DCMAKE_PREFIX_PATH=${cmake_prefix_path}"
                "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
                "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
                "-DCMAKE_INSTALL_PREFIX=${install_dir}"
                "-DCMAKE_BUILD_TYPE=release"
                ${arg_CMAKE_ARGS})

            ExternalProject_Add(${target}
                SOURCE_DIR ${src_dir}
                INSTALL_DIR "${install_dir}"
                STEP_TARGETS build
                DEPENDS ${arg_DEPENDS}
                LIST_SEPARATOR :
                CMAKE_ARGS "${cmake_args}"
            )

            ExternalProject_Add_Step(${target} reconfigure
                COMMAND "${CMAKE_COMMAND}" -E touch "${reconfigure_tag}"
                DEPENDEES install
            )
        else()
            add_error_target(${sub}
                "Building ${sub} library"
                "The git submodule for ${sub} is not available")
        endif()
    endif()
endfunction()

# Catch library for contra

find_package(Catch2 QUIET)
add_external_cmake_submodule(catch
    TARGET make-catch
    UNLESS Catch2_FOUND
    CMAKE_ARGS
    "-DCATCH_BUILD_TESTING=OFF"
    "-DCATCH_BUILD_EXAMPLES=OFF"
    "-DCATCH_INSTALL_DOCS=OFF"
    "-DCATCH_INSTALL_HELPERS=OFF"
)

# Conduit library for contra and nesci
# Note: Findconduit.cmake provided in toplevel cmake subdirectory.

find_package(conduit QUIET)
add_external_cmake_submodule(conduit
    TARGET make-conduit
    UNLESS conduit_FOUND
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/conduit/src"
    CMAKE_ARGS
    "-DENABLE_TESTS=OFF"
    "-DENABLE_DOCS=OFF"
    "-DCMAKE_SKIP_RPATH=ON"
)

# ZeroMQ C++ binding for contra.

find_package(cppzmq QUIET)
add_external_cmake_submodule(cppzmq
    TARGET make-cppzmq
    UNLESS cppzmq_FOUND
    CMAKE_ARGS
    "-DENABLE_DRAFTS=OFF"
    "-DCPPZMQ_BUILD_TESTS=OFF"
)

# Contra library

#find_package(contra QUIET)
find_package(contra)
add_external_cmake_submodule(contra
    TARGET make-contra
    UNLESS contra_FOUND
    DEPENDS make-cppzmq make-conduit make-catch
    CMAKE_ARGS
    "-DENABLE_DRAFTS=OFF"
    "-DCPPZMQ_BUILD_TESTS=OFF"
)

# Nesci library

find_package(nesci COMPONENTS layout producer)
add_external_cmake_submodule(nesci
    TARGET make-nesci
    UNLESS nesci_FOUND
    DEPENDS make-conduit make-catch
)

# Make aliases for contra and nesci for use by examples.

if(contra_FOUND)
    add_library(ext-contra ALIAS contra)
    message(STATUS "*** DEBUG: made contra alias")
else()
    # Running make-contra should trigger a reconfigure.
    add_library(ext-contra INTERFACE)
    add_dependencies(ext-contra make-contra)
endif()

include("${reconfigure_tag}")
if(nesci_FOUND)
    add_library(ext-nesci ALIAS nesci_producer)
    message(STATUS "*** DEBUG: made nesci alias")
else()
    # Running make-nesci should trigger a reconfigure.
    add_library(ext-nesci INTERFACE)
    add_dependencies(ext-nesci make-nesci)
endif()


