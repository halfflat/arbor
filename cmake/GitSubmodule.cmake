# Call to ensure that the git submodule in location `path` is loaded.
# If the submodule is not loaded, an error message that describes
# how to update the submodules is printed.
# Sets the variable given in success_var to `ON` if the submodule is available,
# or `OFF` otherwise.

function(check_git_submodule success_var path)
    set(${success_var} ON PARENT_SCOPE)

    get_filename_component(dotgit "${path}/.git" ABSOLUTE)
    get_filename_component(name "${path}" NAME)
    if(NOT EXISTS ${dotgit} AND NOT "QUIET" IN_LIST ARGN)
        message(
            "\nThe git submodule for ${name} is not available.\n"
            "To check out all submodules use the following commands:\n"
            "    git submodule init\n"
            "    git submodule update\n"
            "Or download submodules recursively when checking out:\n"
            "    git clone --recursive https://github.com/eth-cscs/arbor.git\n"
        )

        # if the repository was not available, and git failed, set AVAIL to false
        set(${success_var} OFF PARENT_SCOPE)
    endif()
endfunction()
