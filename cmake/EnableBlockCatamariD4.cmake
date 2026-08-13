if(NOT DEFINED SOURCE_DIR)
    message(FATAL_ERROR "EnableBlockCatamariD4.cmake requires SOURCE_DIR")
endif()
if(NOT DEFINED COMPONENT)
    message(FATAL_ERROR "EnableBlockCatamariD4.cmake requires COMPONENT")
endif()

function(meshfem_replace_exact file_path old_text new_text description)
    if(NOT EXISTS "${file_path}")
        message(FATAL_ERROR "Cannot apply ${description}: missing ${file_path}")
    endif()

    file(READ "${file_path}" contents)
    string(FIND "${contents}" "${new_text}" new_position)
    if(NOT new_position EQUAL -1)
        message(STATUS "BlockCatamari d=4 patch already present: ${description}")
        return()
    endif()

    string(FIND "${contents}" "${old_text}" old_position)
    if(old_position EQUAL -1)
        message(FATAL_ERROR
            "Cannot apply ${description}; the dependency source has drifted. "
            "Expected text was not found in ${file_path}.")
    endif()

    string(REPLACE "${old_text}" "${new_text}" updated_contents "${contents}")
    if(updated_contents STREQUAL contents)
        message(FATAL_ERROR "Failed to apply ${description} to ${file_path}")
    endif()

    file(WRITE "${file_path}" "${updated_contents}")
    message(STATUS "Applied BlockCatamari d=4 patch: ${description}")
endfunction()

if(COMPONENT STREQUAL "MeshFEMSparse")
    meshfem_replace_exact(
        "${SOURCE_DIR}/src/lib/MeshFEMSparse/Solvers/CatamariConverter.hh"
        [=[#define MAX_INSTANTIATED_BLOCK_SIZE 3]=]
        [=[#define MAX_INSTANTIATED_BLOCK_SIZE 4]=]
        "raise MeshFEMSparse's instantiated block-size limit")

    meshfem_replace_exact(
        "${SOURCE_DIR}/src/lib/MeshFEMSparse/Solvers/CatamariFactorizer.hh"
        [=[    void setUseBlockAccel(bool u) { m_useBlockAccel = u; }
    bool getUseBlockAccel() const { return m_useBlockAccel; }

    // For benchmarking comparisons: disable the use of block accelerations only for numeric factorization or solves.]=]
        [=[    void setUseBlockAccel(bool u) { m_useBlockAccel = u; }
    bool getUseBlockAccel() const { return m_useBlockAccel; }

    // Actual uniform block size selected during symbolic factorization.
    // A value of 1 means the factorizer is using the scalar fallback.
    size_t getFactorizationBlockSize() const { return m_blockSize; }

    // For benchmarking comparisons: disable the use of block accelerations only for numeric factorization or solves.]=]
        "expose the selected factorization block size for regression testing")
elseif(COMPONENT STREQUAL "BlockCatamari")
    meshfem_replace_exact(
        "${SOURCE_DIR}/include/catamari/sparse_ldl/supernodal/factorization.hpp"
        [=[          if      (blockSize == 3) result = BlockLeftLooking<3>();
          else if (blockSize == 2) result = BlockLeftLooking<2>();
          else                     result = BlockLeftLooking<1>();]=]
        [=[          if      (blockSize == 4) result = BlockLeftLooking<4>();
          else if (blockSize == 3) result = BlockLeftLooking<3>();
          else if (blockSize == 2) result = BlockLeftLooking<2>();
          else                     result = BlockLeftLooking<1>();]=]
        "dispatch d=4 left-looking numeric factorization")

    meshfem_replace_exact(
        "${SOURCE_DIR}/include/catamari/sparse_ldl/supernodal/factorization.hpp"
        [=[      if (blockSize == 3) return BlockRightLooking<3>();
      if (blockSize == 2) return BlockRightLooking<2>();
      return BlockRightLooking<1>();]=]
        [=[      if (blockSize == 4) return BlockRightLooking<4>();
      if (blockSize == 3) return BlockRightLooking<3>();
      if (blockSize == 2) return BlockRightLooking<2>();
      return BlockRightLooking<1>();]=]
        "dispatch d=4 right-looking numeric factorization")

    meshfem_replace_exact(
        "${SOURCE_DIR}/include/catamari/sparse_ldl/supernodal/factorization/solve-impl.hpp"
        [=[        if (block_size == 3) {
            OpenMPLowerTriangularSolve<3>(&permuted_right_hand_sides, &shared_state);
            OpenMPDiagonalSolve(&permuted_right_hand_sides);
            OpenMPLowerTransposeTriangularSolve<3>(&permuted_right_hand_sides, &shared_state);
        }
        else if (block_size == 2) {]=]
        [=[        if (block_size == 4) {
            OpenMPLowerTriangularSolve<4>(&permuted_right_hand_sides, &shared_state);
            OpenMPDiagonalSolve(&permuted_right_hand_sides);
            OpenMPLowerTransposeTriangularSolve<4>(&permuted_right_hand_sides, &shared_state);
        }
        else if (block_size == 3) {
            OpenMPLowerTriangularSolve<3>(&permuted_right_hand_sides, &shared_state);
            OpenMPDiagonalSolve(&permuted_right_hand_sides);
            OpenMPLowerTransposeTriangularSolve<3>(&permuted_right_hand_sides, &shared_state);
        }
        else if (block_size == 2) {]=]
        "dispatch d=4 parallel triangular solves")

    meshfem_replace_exact(
        "${SOURCE_DIR}/include/catamari/sparse_ldl/supernodal/factorization/solve-impl.hpp"
        [=[      if (block_size == 3) {
          LowerTriangularSolve<3>(&permuted_right_hand_sides);
          DiagonalSolve(&permuted_right_hand_sides);
          LowerTransposeTriangularSolve<3>(&permuted_right_hand_sides);
      }
      else if (block_size == 2) {]=]
        [=[      if (block_size == 4) {
          LowerTriangularSolve<4>(&permuted_right_hand_sides);
          DiagonalSolve(&permuted_right_hand_sides);
          LowerTransposeTriangularSolve<4>(&permuted_right_hand_sides);
      }
      else if (block_size == 3) {
          LowerTriangularSolve<3>(&permuted_right_hand_sides);
          DiagonalSolve(&permuted_right_hand_sides);
          LowerTransposeTriangularSolve<3>(&permuted_right_hand_sides);
      }
      else if (block_size == 2) {]=]
        "dispatch d=4 serial triangular solves")
else()
    message(FATAL_ERROR "Unknown COMPONENT '${COMPONENT}'")
endif()
