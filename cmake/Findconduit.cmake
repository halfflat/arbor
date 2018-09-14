# Hack up a version of Findconduit that publishes include directory too.

if(NOT conduit_FOUND)
    find_file(conduit_cmake PATH_SUFFIXES lib/cmake share/cmake NAMES conduit.cmake)
    if(conduit_cmake)
        include(${conduit_cmake})
        if(TARGET conduit)
            get_target_property(_location conduit LOCATION)
            get_filename_component(_location "${_location}" PATH)
            get_filename_component(_location "${_location}" PATH)

            find_file(_header HINTS "${location_}" NAMES conduit/conduit.hpp conduit/conduit.h)
            if(_header)
                get_filename_component(_header "${_header}" PATH)
                get_filename_component(_header "${_header}" PATH)
            endif()

            foreach(conduit_target conduit conduit_relay conduit_blueprint)
                if(TARGET ${conduit_target})
                    set_target_properties(${conduit_target} PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${_header}")
                endif()
            endforeach()
        endif()
    endif()

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(conduit DEFAULT_MSG conduit_cmake)
endif()
