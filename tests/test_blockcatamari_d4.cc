#include <MeshFEMCore/Parallelism.hh>
#include <MeshFEMSparse/SystemAssembler.hh>
#include <MeshFEMSparse/Solvers/CatamariFactorizer.hh>

#include <array>
#include <vector>

// WARNING: catch2/catch.hpp defines a BENCHMARK macro, so include it after
// MeshFEM headers.
#include <catch2/catch.hpp>

using namespace MeshFEM;

namespace {
struct TBBThreadLimit {
    explicit TBBThreadLimit(int numThreads) {
        set_max_num_tbb_threads(numThreads);
    }
    ~TBBThreadLimit() {
        unset_max_num_tbb_threads();
    }
};

struct MatrixCase {
    const char *name;
    const BlockCSCHessianBase *matrix;
};

struct PinCase {
    const char *name;
    const std::vector<size_t> *fixedVars;
};
}

TEST_CASE("BlockCatamari factors and solves uniform 4D blocks",
          "[blockcatamari_d4]") {
    constexpr size_t blockSize = 4;
    constexpr size_t numBlocks = 4;
    constexpr double shift = 0.25;
    constexpr double transportWeight = 0.2;

    const std::array<std::array<size_t, 2>, numBlocks - 1> stencils {{
        {{0, 1}},
        {{1, 2}},
        {{2, 3}},
    }};

    SystemAssembler<blockSize> assembler(numBlocks);
    auto H = assembler.blockSparsityPattern(
        stencils.size(),
        [&stencils](size_t ei) { return stencils[ei]; });

    // Start with I and add local PSD transport terms
    //   w ||x_{i + 1} - R x_i||^2.
    // The dense R deliberately couples all four coordinates, exercising
    // actual 4x4 block arithmetic instead of four decoupled scalar systems.
    H->setIdentity(/* preserveSparsity = */ true);
    Eigen::Matrix4d transport;
    transport <<
         1.00,  0.08, -0.03,  0.02,
        -0.04,  0.97,  0.06, -0.01,
         0.05, -0.02,  1.03,  0.07,
         0.01,  0.04, -0.05,  0.99;
    const Eigen::Matrix4d sourceContribution =
        transportWeight * transport.transpose() * transport;

    for (size_t bi = 0; bi + 1 < numBlocks; ++bi) {
        const size_t source = blockSize * bi;
        const size_t target = blockSize * (bi + 1);
        for (size_t a = 0; a < blockSize; ++a) {
            for (size_t b = a; b < blockSize; ++b) {
                H->addNZScalar(source + a, source + b,
                               sourceContribution(a, b));
            }
            H->addNZScalar(target + a, target + a, transportWeight);
            for (size_t b = 0; b < blockSize; ++b) {
                H->addNZScalar(source + a, target + b,
                               -transportWeight * transport(b, a));
            }
        }
    }

    const auto scalarH = H->toScalar();
    auto compressedH = BlockCSCHessianFromScalar(scalarH, blockSize);
    const std::array<MatrixCase, 2> matrices {{
        {"assembled", H.get()},
        {"scalar-compressed", compressedH.get()},
    }};

    const std::vector<size_t> noFixedVars;
    const std::vector<size_t> fixedFirstBlock {0, 1, 2, 3};
    const std::array<PinCase, 2> pinCases {{
        {"none", &noFixedVars},
        {"whole-d4-block", &fixedFirstBlock},
    }};

    Eigen::VectorXd expected = Eigen::VectorXd::LinSpaced(
        blockSize * numBlocks, 1.0, double(blockSize * numBlocks));
    expected.head(blockSize).setZero();
    const Eigen::VectorXd rhs = H->apply(expected);
    const Eigen::VectorXd shiftedRhs = rhs + shift * expected;

    for (const MatrixCase &matrixCase : matrices) {
        for (const PinCase &pinCase : pinCases) {
            // Cover both factorization implementations and solve thread budgets.
            for (bool useLeftLooking : {false, true}) {
                for (int numThreads : {1, 2}) {
                    TBBThreadLimit threadLimit(numThreads);

                    CatamariFactorizer factorizer;
                    factorizer.orderingMethod =
                        CatamariFactorizer::OrderingMethod::AMD;
                    factorizer.setUseLeftLooking(useLeftLooking);
                    factorizer.factorizeSymbolic(
                        *matrixCase.matrix, *pinCase.fixedVars);

                    INFO("matrix = " << matrixCase.name);
                    INFO("pins = " << pinCase.name);
                    INFO("left-looking = " << useLeftLooking);
                    INFO("threads = " << numThreads);
                    REQUIRE(factorizer.getFactorizationBlockSize() == blockSize);

                    factorizer.factorizeNumeric(*matrixCase.matrix);
                    REQUIRE((factorizer.solve(rhs) - expected).norm()
                            / expected.norm() < 1e-11);

                    // Exercise MORSE's Levenberg--Marquardt A + tau I path while
                    // retaining the same symbolic factorization and d=4 block layout.
                    factorizer.factorizeNumericWithShift(
                        *matrixCase.matrix, shift);
                    REQUIRE((factorizer.solve(shiftedRhs) - expected).norm()
                            / expected.norm() < 1e-11);
                }
            }
        }
    }
}
