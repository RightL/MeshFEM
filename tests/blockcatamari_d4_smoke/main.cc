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

struct MatrixCase {
    const char *name;
    const BlockCSCHessianBase *matrix;
};

struct PinCase {
    const char *name;
    const std::vector<size_t> *fixedVars;
};

bool checkSolution(const Eigen::VectorXd &actual,
                   const Eigen::VectorXd &expected,
                   const char *matrixKind,
                   const char *pinKind,
                   bool useLeftLooking,
                   int numThreads,
                   const char *solveKind) {
    const double relativeError =
        (actual - expected).norm() / expected.norm();
    if (relativeError < 1e-11) return true;

    std::cerr << solveKind << " solve failed for matrix=" << matrixKind
              << ", pins=" << pinKind
              << ", left-looking=" << useLeftLooking
              << ", threads=" << numThreads
              << ", relative error=" << relativeError << '\n';
    return false;
}
}

int main() {
    try {
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

        // Exercise both a natively assembled d=4 matrix and the migration path
        // from an existing scalar SuiteSparseMatrix representation.
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
                for (bool useLeftLooking : {false, true}) {
                    for (int numThreads : {1, 2}) {
                        TBBThreadLimit threadLimit(numThreads);

                        CatamariFactorizer factorizer;
                        factorizer.orderingMethod =
                            CatamariFactorizer::OrderingMethod::AMD;
                        factorizer.setUseLeftLooking(useLeftLooking);
                        factorizer.factorizeSymbolic(
                            *matrixCase.matrix, *pinCase.fixedVars);

                        if (factorizer.getFactorizationBlockSize() != blockSize) {
                            std::cerr << "Expected d=4 block factorization for "
                                      << matrixCase.name << " with pins="
                                      << pinCase.name << ", got d="
                                      << factorizer.getFactorizationBlockSize()
                                      << '\n';
                            return 1;
                        }

                        factorizer.factorizeNumeric(*matrixCase.matrix);
                        if (!checkSolution(factorizer.solve(rhs), expected,
                                           matrixCase.name, pinCase.name,
                                           useLeftLooking, numThreads,
                                           "Unshifted")) {
                            return 1;
                        }

                        // This is the Levenberg--Marquardt path used by MORSE:
                        // refactor the same pattern after adding tau I.
                        factorizer.factorizeNumericWithShift(
                            *matrixCase.matrix, shift);
                        if (!checkSolution(factorizer.solve(shiftedRhs), expected,
                                           matrixCase.name, pinCase.name,
                                           useLeftLooking, numThreads,
                                           "Shifted")) {
                            return 1;
                        }
                    }
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
