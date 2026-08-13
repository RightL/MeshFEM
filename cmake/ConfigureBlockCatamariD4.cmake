function(meshfem_configure_blockcatamari_d4)
    set(patch_script "${MESHFEM_ROOT}/cmake/EnableBlockCatamariD4.cmake")

    # FetchContent keeps the first declaration it sees. Declaring Catamari here
    # lets MeshFEMSparse's later nested FetchContent_MakeAvailable() use the
    # patched source without requiring a separate fork of the dependency.
    set(catamari_source_dir "${MESHFEM_EXTERNAL}/catamari")
    if(EXISTS "${catamari_source_dir}/CMakeLists.txt")
        execute_process(
            COMMAND ${CMAKE_COMMAND}
                "-DSOURCE_DIR=${catamari_source_dir}"
                "-DCOMPONENT=BlockCatamari"
                -P "${patch_script}"
            COMMAND_ERROR_IS_FATAL ANY
        )
    endif()

    FetchContent_Declare(catamari
        GIT_REPOSITORY https://github.com/MeshFEM/BlockCatamari.git
        GIT_TAG        master
        SOURCE_DIR     ${catamari_source_dir}
        PATCH_COMMAND
            ${CMAKE_COMMAND}
            "-DSOURCE_DIR=<SOURCE_DIR>"
            "-DCOMPONENT=BlockCatamari"
            -P "${patch_script}"
    )

    # Apply the same idempotent patch immediately for an existing editable
    # checkout. For a fresh checkout, FetchContent runs PATCH_COMMAND after
    # cloning and before adding the dependency as a subdirectory.
    set(meshfemsparse_source_dir "${MESHFEM_EXTERNAL}/MeshFEMSparse")
    if(EXISTS "${meshfemsparse_source_dir}/CMakeLists.txt")
        execute_process(
            COMMAND ${CMAKE_COMMAND}
                "-DSOURCE_DIR=${meshfemsparse_source_dir}"
                "-DCOMPONENT=MeshFEMSparse"
                -P "${patch_script}"
            COMMAND_ERROR_IS_FATAL ANY
        )
    endif()

    FetchContent_Declare(MeshFEMSparse
        GIT_REPOSITORY ${MESHFEMSPARSE_GIT_REPOSITORY}
        GIT_TAG        ${MESHFEMSPARSE_GIT_TAG}
        SOURCE_DIR     ${meshfemsparse_source_dir}
        PATCH_COMMAND
            ${CMAKE_COMMAND}
            "-DSOURCE_DIR=<SOURCE_DIR>"
            "-DCOMPONENT=MeshFEMSparse"
            -P "${patch_script}"
    )
endfunction()
