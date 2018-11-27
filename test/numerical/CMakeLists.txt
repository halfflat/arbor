set(numerical_mechs
    decay_injection
)

set(external_modcc)
if(ARB_WITH_EXTERNAL_MODCC)
    set(external_modcc MODCC ${modcc})
endif()
set(mech_dir ${CMAKE_CURRENT_BINARY_DIR}/mechanisms)

build_modules(
    ${numerical_mechs}
    SOURCE_DIR mod
    DEST_DIR "${mech_dir}"
    ${external_modcc}
    MODCC_FLAGS -t cpu ${ARB_MODCC_FLAGS}
    GENERATES .hpp _cpu.cpp
    TARGET numerical_mechs
)

set(numerical_mech_sources)
foreach(mech ${numerical_mechs})
    list(APPEND numerical_mech_sources ${mech_dir}/${mech}_cpu.cpp)
endforeach()

add_executable(rc-circuit EXCLUDE_FROM_ALL rc-circuit.cpp ${numerical_mech_sources})
target_link_libraries(rc-circuit PRIVATE arbor arbor-sup)
add_dependencies(rc-circuit numerical_mechs)
