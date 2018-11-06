add_executable(rc-circuit EXCLUDE_FROM_ALL rc-circuit.cpp)

target_link_libraries(rc-circuit PRIVATE arbor arbor-sup)
