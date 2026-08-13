cmake_minimum_required(VERSION 3.14...3.31)
project(BlockCatamariD4Smoke LANGUAGES CXX)

include(CTest)
include(FetchContent)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

get_filename_component(MESHFEM_ROOT "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)
set(MESHFEM_EXTERNAL "${CMAKE_BINARY_DIR}/external")
list(APPEND CMAKE_MODULE_PATH "${MESHFEM_ROOT}/cmake")

set(MESHFEM_NATIVE OFF CACHE BOOL "" FORCE)
set(MESHFEMCORE_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(MESHFEMSPARSE_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(MESHFEM_WITH_CATAMARI ON CACHE BOOL "" FORCE)
set(MESHFEM_USE_LEGACY_CATAMARI OFF CACHE BOOL "" FORCE)
set(MESHFEM_WITH_CHOLMOD ON CACHE BOOL "" FORCE)
set(MESHFEM_WITH_UMFPACK OFF CACHE BOOL "" FORCE)
set(MESHFEM_WITH_SCOTCH OFF CACHE BOOL "" FORCE)

# BlockCatamari respects an existing BLAS::BLAS target. This avoids its MKL-
# specific default search on generic CI machines. Its dense Cholesky kernels
# also call LAPACK directly; generic BLAS packages do not bundle those symbols.
set(BLA_VENDOR Generic)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

set(MESHFEMCORE_GIT_REPOSITORY "https://github.com/MeshFEM/MeshFEMCore.git")
set(MESHFEMCORE_GIT_TAG "main")
FetchContent_Declare(MeshFEMCore
    GIT_REPOSITORY ${MESHFEMCORE_GIT_REPOSITORY}
    GIT_TAG        ${MESHFEMCORE_GIT_TAG}
    SOURCE_DIR     ${MESHFEM_EXTERNAL}/MeshFEMCore
)
FetchContent_MakeAvailable(MeshFEMCore)
FetchContent_GetProperties(MeshFEMCore)
list(APPEND CMAKE_MODULE_PATH "${meshfemcore_SOURCE_DIR}/cmake")
list(REMOVE_DUPLICATES CMAKE_MODULE_PATH)

set(MESHFEMSPARSE_GIT_REPOSITORY "https://github.com/MeshFEM/MeshFEMSparse.git")
set(MESHFEMSPARSE_GIT_TAG "main")
include(ConfigureBlockCatamariD4)
meshfem_configure_blockcatamari_d4()
FetchContent_MakeAvailable(MeshFEMSparse)

add_executable(blockcatamari_d4_smoke main.cc)
target_link_libraries(blockcatamari_d4_smoke PRIVATE MeshFEMSparse LAPACK::LAPACK)

add_test(NAME blockcatamari_d4_smoke COMMAND blockcatamari_d4_smoke)
