#include <MeshFEMCore/Parallelism.hh>
#include <MeshFEMSparse/SystemAssembler.hh>
#include <MeshFEMSparse/Solvers/CatamariFactorizer.hh>

#include <array>
#include <vector>

// WARNING: catch2/catch.hpp defines a BENCHMARK macro, so include it after
// MeshFEM headers.
#include <catch2/catch.hpp>

using namespace MeshFEM;

TEST_CASE("BlockCatamari factors and solves uniform 4D blocks",
          "[blockcatamari_d4]") {
    constexpr size_t blockSize = 4;
    constexpr size_t numBlocks = 4;

    const std::array<std::array<size_t, 2>, numBlocks - 1> stencils {{
        {{0, 1}},
        {{1, 2}},
        {{2, 3}},
    }};

    SystemAssembler<blockSize> assembler(numBlocks);
    auto H = assembler.blockSparsityPattern(
        stencils.size(),
        [&stencils](size_t ei) { return stencils[ei]; });

    // Block-tridiagonal SPD matrix with a nontrivial 4x4 block pattern.
    H->setIdentity(/* preserveSparsity = */ true);
    for (size_t bi = 0; bi + 1 < numBlocks; ++bi) {
        for (size_t c = 0; c < blockSize; ++c) {
            H->addNZScalar(blockSize * bi + c,
                           blockSize * (bi + 1) + c,
                           -0.1);
        }
    }

    CatamariFactorizer factorizer;
    factorizer.orderingMethod = CatamariFactorizer::OrderingMethod::AMD;

    // MORSE fixes all four coordinates of a boundary field value together.
    const std::vector<size_t> fixedVars {0, 1, 2, 3};
    factorizer.factorizeSymbolic(*H, fixedVars);
    REQUIRE(factorizer.getFactorizationBlockSize() == blockSize);

    factorizer.factorizeNumeric(*H);

    Eigen::VectorXd expected = Eigen::VectorXd::LinSpaced(
        blockSize * numBlocks, 1.0, double(blockSize * numBlocks));
    expected.head(blockSize).setZero();

    const Eigen::VectorXd rhs = H->apply(expected);
    const Eigen::VectorXd actual = factorizer.solve(rhs);

    REQUIRE((actual - expected).norm() / expected.norm() < 1e-11);
}
