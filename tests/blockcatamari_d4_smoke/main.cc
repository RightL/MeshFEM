#include <MeshFEMCore/Parallelism.hh>
#include <MeshFEMSparse/SystemAssembler.hh>
#include <MeshFEMSparse/Solvers/CatamariFactorizer.hh>

#include <array>
#include <exception>
#include <iostream>
#include <vector>

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

bool checkSolution(const Eigen::VectorXd &actual,
                   const Eigen::VectorXd &expected,
                   bool useLeftLooking,
                   int numThreads,
                   const char *solveKind) {
    const double relativeError =
        (actual - expected).norm() / expected.norm();
    if (relativeError < 1e-11) return true;

    std::cerr << solveKind << " solve failed for left-looking="
              << useLeftLooking << ", threads=" << numThreads
              << ", relative error=" << relativeError << '\n';
    return false;
}
}

int main() {
    try {
        constexpr size_t blockSize = 4;
        constexpr size_t numBlocks = 4;
        constexpr double shift = 0.25;

        const std::array<std::array<size_t, 2>, numBlocks - 1> stencils {{
            {{0, 1}},
            {{1, 2}},
            {{2, 3}},
        }};

        SystemAssembler<blockSize> assembler(numBlocks);
        auto H = assembler.blockSparsityPattern(
            stencils.size(),
            [&stencils](size_t ei) { return stencils[ei]; });

        H->setIdentity(/* preserveSparsity = */ true);
        for (size_t bi = 0; bi + 1 < numBlocks; ++bi) {
            for (size_t c = 0; c < blockSize; ++c) {
                H->addNZScalar(blockSize * bi + c,
                               blockSize * (bi + 1) + c,
                               -0.1);
            }
        }

        Eigen::VectorXd expected = Eigen::VectorXd::LinSpaced(
            blockSize * numBlocks, 1.0, double(blockSize * numBlocks));
        expected.head(blockSize).setZero();
        const Eigen::VectorXd rhs = H->apply(expected);
        const Eigen::VectorXd shiftedRhs = rhs + shift * expected;
        const std::vector<size_t> fixedVars {0, 1, 2, 3};

        for (bool useLeftLooking : {false, true}) {
            for (int numThreads : {1, 2}) {
                TBBThreadLimit threadLimit(numThreads);

                CatamariFactorizer factorizer;
                factorizer.orderingMethod = CatamariFactorizer::OrderingMethod::AMD;
                factorizer.setUseLeftLooking(useLeftLooking);
                factorizer.factorizeSymbolic(*H, fixedVars);

                if (factorizer.getFactorizationBlockSize() != blockSize) {
                    std::cerr << "Expected d=4 block factorization, got d="
                              << factorizer.getFactorizationBlockSize() << '\n';
                    return 1;
                }

                factorizer.factorizeNumeric(*H);
                if (!checkSolution(factorizer.solve(rhs), expected,
                                   useLeftLooking, numThreads, "Unshifted")) {
                    return 1;
                }

                // This is the Levenberg--Marquardt path used by MORSE:
                // refactor the same pattern after adding tau I.
                factorizer.factorizeNumericWithShift(*H, shift);
                if (!checkSolution(factorizer.solve(shiftedRhs), expected,
                                   useLeftLooking, numThreads, "Shifted")) {
                    return 1;
                }
            }
        }

        std::cout << "BlockCatamari d=4 smoke test passed\n";
        return 0;
    }
    catch (const std::exception &e) {
        std::cerr << "BlockCatamari d=4 smoke test threw: " << e.what() << '\n';
        return 1;
    }
}
