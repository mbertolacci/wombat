#include <RcppEigen.h>

using Rcpp::IntegerVector;
using Rcpp::List;
using Rcpp::NumericVector;
using Rcpp::S4;

using Eigen::MatrixXd;
typedef Eigen::SparseMatrix<double> MatrixSd;
typedef Eigen::Map<MatrixSd> MapSd;
typedef Eigen::Map<MatrixXd> MapXd;

void checkCompatibility(const S4& a, const S4& b) {
    IntegerVector aDim = a.slot("Dim");
    IntegerVector bDim = b.slot("Dim");

    if (!Rcpp::is_true(Rcpp::all(aDim == bDim))) {
        Rcpp::stop("Matrices have incompatible dimensions");
    }
}

bool isMatrixSparse(const S4& a) {
    return a.is("CsparseMatrix");
}

bool isMatrixDense(const S4& a) {
    return a.is("dgeMatrix") || a.is("dsyMatrix");
}

bool isMatrixDiagonal(const S4& a) {
    return a.is("ddiMatrix");
}

bool isMatrixSymmetric(const S4& a) {
    return a.is("dsyMatrix") || a.is("dsCMatrix");
}

bool isMatrixUp(const S4& a) {
    if (isMatrixSymmetric(a)) {
        std::string uplo = a.slot("uplo");
        return uplo == "U";
    }
    return false;
}

MapSd toEigenSparse(S4 x) {
    IntegerVector xDimensions = x.slot("Dim");

    IntegerVector xP = x.slot("p");
    IntegerVector xI = x.slot("i");
    NumericVector xX = x.slot("x");

    return MapSd(
        xDimensions[0],
        xDimensions[1],
        xP[xDimensions[1]],
        xP.begin(),
        xI.begin(),
        xX.begin()
    );
}

MapXd toEigenDense(S4 x) {
    NumericVector xX = x.slot("x");
    Rcpp::Dimension xDim = x.slot("Dim");
    return MapXd(
        xX.begin(),
        xDim[0],
        xDim[1]
    );
}

S4 toMatrixSparse(MatrixSd& x, bool symmetric, bool up) {
    x.makeCompressed();

    S4 output(symmetric ? "dsCMatrix" : "dgCMatrix");
    const int numNonZero = x.nonZeros();
    output.slot("Dim") = Rcpp::Dimension(x.rows(), x.cols());
    output.slot("i") = IntegerVector(x.innerIndexPtr(), x.innerIndexPtr() + numNonZero);
    output.slot("p") = IntegerVector(x.outerIndexPtr(), x.outerIndexPtr() + x.outerSize() + 1);
    output.slot("x") = NumericVector(x.valuePtr(), x.valuePtr() + numNonZero);
    if (symmetric) {
        output.slot("uplo") = up ? "U" : "L";
    }
    return output;
}

S4 toMatrixDense(const MatrixXd& x, bool symmetric, bool up) {
    S4 output(symmetric ? "dsyMatrix" : "dgeMatrix");
    output.slot("x") = NumericVector(x.data(), x.data() + x.size());
    output.slot("Dim") = IntegerVector({
        static_cast<int>(x.rows()),
        static_cast<int>(x.cols())
    });
    output.slot("Dimnames") = List({ R_NilValue, R_NilValue });
    output.slot("factors") = List();
    if (symmetric) {
        output.slot("uplo") = up ? "U" : "L";
    }
    return output;
}

template <typename Derived>
typename Derived::PlainObject sparseSymmetricToGeneral(Eigen::SparseMatrixBase<Derived>& x, bool up) {
    typedef typename Derived::PlainObject Output;

    Output output(x.rows(), x.cols());
    if (up) {
        output = x.template selfadjointView<Eigen::Upper>();
    } else {
        output = x.template selfadjointView<Eigen::Lower>();
    }

    return output;
}

S4 fastAddSparseSparse(const S4& aR, const S4& bR) {
    bool symmetric = false;
    bool transposeB = false;
    if (isMatrixSymmetric(aR) && isMatrixSymmetric(bR)) {
        symmetric = true;
        std::string aUplo = aR.slot("uplo");
        std::string bUplo = bR.slot("uplo");
        if (aUplo != bUplo) {
            transposeB = true;
        }
    }
    MapSd a = toEigenSparse(aR);
    MapSd b = toEigenSparse(bR);

    MatrixSd result;
    if (transposeB) {
        result = a + MatrixSd(b.transpose());
    } else {
        if (!symmetric && isMatrixSymmetric(aR)) {
            result = sparseSymmetricToGeneral(a, isMatrixUp(aR)) + b;
        } else if (!symmetric && isMatrixSymmetric(bR)) {
            result = a + sparseSymmetricToGeneral(b, isMatrixUp(bR));
        } else {
            result = a + b;
        }
    }

    return toMatrixSparse(
        result,
        symmetric,
        isMatrixUp(aR)
    );
}

S4 fastAddSparseDiagonal(const S4& aR, const S4& bR) {
    bool symmetric = false;
    bool up = false;
    if (isMatrixSymmetric(aR)) {
        symmetric = true;
        std::string aUplo = aR.slot("uplo");
        if (aUplo == "U") up = true;
    }

    MatrixSd result = toEigenSparse(aR);

    std::string bDiag = bR.slot("diag");
    if (bDiag == "U") {
        MatrixSd b(result.rows(), result.rows());
        b.setIdentity();

        result += b;
    } else {
        NumericVector bX = bR.slot("x");
        Eigen::Map<Eigen::VectorXd> b = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(bX);

        result += MatrixSd(b.asDiagonal());
    }
    return toMatrixSparse(result, symmetric, up);
}

S4 fastAddSparseDense(const S4& aR, const S4& bR) {
    MapSd a = toEigenSparse(aR);
    MapXd b = toEigenDense(bR);

    MatrixXd result;
    if (isMatrixSymmetric(aR) && isMatrixSymmetric(bR)) {
        if (isMatrixUp(aR) && isMatrixUp(bR)) {
            result = b;
        } else {
            result = b.transpose();
        }
        result += a;
    } else {
        if (isMatrixSymmetric(bR)) {
            if (isMatrixUp(bR)) {
                result = b.selfadjointView<Eigen::Upper>();
            } else {
                result = b.selfadjointView<Eigen::Lower>();
            }
        } else {
            result = b;
        }

        if (isMatrixSymmetric(aR)) {
            result += sparseSymmetricToGeneral(a, isMatrixUp(aR));
        } else {
            result += a;
        }
    }

    return toMatrixDense(
        result,
        isMatrixSymmetric(aR) && isMatrixSymmetric(bR),
        isMatrixUp(aR)
    );
}

S4 fastAddDiagonalDiagonal(const S4& aR, const S4& bR) {
    IntegerVector aDim = aR.slot("Dim");
    int n = aDim[0];
    std::string aDiag = aR.slot("diag");
    std::string bDiag = bR.slot("diag");

    NumericVector aX = aR.slot("x");
    NumericVector bX = bR.slot("x");
    if (aDiag == "U") {
        aX = NumericVector(n, 1.0);
    }
    if (bDiag == "U") {
        bX = NumericVector(n, 1.0);
    }

    S4 result = S4("ddiMatrix");
    result.slot("Dim") = aDim;
    result.slot("Dimnames") = aR.slot("Dimnames");
    result.slot("x") = aX + bX;
    result.slot("diag") = "N";
    return result;
}

S4 fastAddDiagonalDense(const S4& aR, const S4& bR) {
    MatrixXd result = toEigenDense(bR);

    std::string aDiag = aR.slot("diag");
    if (aDiag == "U") {
        result.diagonal().array() += 1;
    } else {
        NumericVector aX = aR.slot("x");
        for (int i = 0; i < result.rows(); ++i) {
            result(i, i) += aX[i];
        }
    }
    return toMatrixDense(
        result,
        isMatrixSymmetric(bR),
        isMatrixUp(bR)
    );
}

S4 fastAddDenseDense(const S4& aR, const S4& bR) {
    bool symmetric = isMatrixSymmetric(aR) && isMatrixSymmetric(bR);

    NumericVector aX = aR.slot("x");
    NumericVector bX = bR.slot("x");

    S4 output(symmetric ? "dsyMatrix" : "dgeMatrix");
    output.slot("x") = aX + bX;
    output.slot("Dim") = aR.slot("Dim");
    output.slot("Dimnames") = aR.slot("Dimnames");
    output.slot("factors") = List();
    if (symmetric) {
        output.slot("uplo") = isMatrixUp(aR) ? "U" : "L";
    }

    return output;
}

// Fast addition for general + general and symmetric + symmetric
// [[Rcpp::export(name=".fast_add", rng=false)]]
S4 fast_add(S4 a, S4 b) {
    if (isMatrixDiagonal(a) || isMatrixDense(a)) {
        std::swap(a ,b);
    }

    checkCompatibility(a, b);

    if (isMatrixSparse(a)) {
        if (isMatrixSparse(b)) {
            return fastAddSparseSparse(a, b);
        } else if (isMatrixDiagonal(b)) {
            return fastAddSparseDiagonal(a, b);
        } else if (isMatrixDense(b)) {
            return fastAddSparseDense(a, b);
        }
    } else if (isMatrixDiagonal(a)) {
        if (isMatrixDiagonal(b)) {
            return fastAddDiagonalDiagonal(a, b);
        } else if (isMatrixDense(b)) {
            return fastAddDiagonalDense(a, b);
        }
    } else if (isMatrixDense(a)) {
        if (isMatrixDense(b)) {
            return fastAddDenseDense(a, b);
        } else if (isMatrixDiagonal(b)) {
            return fastAddDiagonalDense(b, a);
        }
    }

    Rcpp::stop("Invalid combination of matrices encountered");
}

// [[Rcpp::export(name=".fast_as_dgeMatrix", rng=false)]]
S4 fast_as_dgeMatrix(const Eigen::MatrixXd& x) {
    return toMatrixDense(x, false, false);
}


// Fast addition for general + general and symmetric + symmetric
// [[Rcpp::export(name=".fast_add2", rng=false)]]
S4 fast_add2(S4 a, S4 b) {
    if (a.inherits("dsCMatrix") && b.inherits("dsCMatrix")) {
        std::string aUplo = a.slot("uplo");
        std::string bUplo = b.slot("uplo");
        if (aUplo != bUplo) {
            Rcpp::stop("a and b must use same upper triangle representation");
        }
    } else if (a.inherits("dsCMatrix") || b.inherits("dsCMatrix")) {
        Rcpp::stop("Symmetric plus general not supported");
    }

    IntegerVector aDims(a.slot("Dim"));
    IntegerVector bDims(b.slot("Dim"));

    if (aDims[0] != bDims[0] || aDims[1] != bDims[1]) {
        Rcpp::stop("Matrix dimensions do not match");
    }

    IntegerVector aI_(a.slot("i"));
    int *aI = &aI_[0];
    IntegerVector aP_(a.slot("p"));
    int *aP = &aP_[0];
    NumericVector aX_(a.slot("x"));
    double *aX = &aX_[0];

    IntegerVector bI_(b.slot("i"));
    int *bI = &bI_[0];
    IntegerVector bP_(b.slot("p"));
    int *bP = &bP_[0];
    NumericVector bX_(b.slot("x"));
    double *bX = &bX_[0];

    int nRows = aDims[0];
    int nColumns = aDims[1];

    // Compute P
    IntegerVector outputP_ = Rcpp::no_init(nColumns + 1);
    int *outputP = &outputP_[0];
    outputP[0] = 0;
    for (int col = 0; col < nColumns; ++col) {
        int aK = aP[col];
        int bK = bP[col];
        outputP[col + 1] = outputP[col];

        while (aK < aP[col + 1] || bK < bP[col + 1]) {
            if (aK == aP[col + 1]) {
                // a has reached end of column, advance b
                ++bK;
            } else if (bK == bP[col + 1]) {
                // b has reached end of column, advance a
                ++aK;
            } else if (aI[aK] < bI[bK]) {
                // a is behind b, advance a
                ++aK;
            } else if (bI[bK] < aI[aK]) {
                // b is behind a, advance b
                ++bK;
            } else {
                // they are at the same row, advance both
                ++aK;
                ++bK;
            }
            ++outputP[col + 1];
        }
    }

    // Compute sum
    IntegerVector outputI_ = Rcpp::no_init(outputP[nColumns]);
    int *outputI = &outputI_[0];
    NumericVector outputX_ = Rcpp::no_init(outputP[nColumns]);
    double *outputX = &outputX_[0];
    for (int col = 0; col < nColumns; ++col) {
        int aK = aP[col];
        int bK = bP[col];

        int outputK = outputP[col];
        while (aK < aP[col + 1] || bK < bP[col + 1]) {
            if (aK == aP[col + 1]) {
                // a has reached end of column, advance b
                outputX[outputK] = bX[bK];
                outputI[outputK] = bI[bK];
                ++bK;
            } else if (bK == bP[col + 1]) {
                // b has reached end of column, advance a
                outputX[outputK] = aX[aK];
                outputI[outputK] = aI[aK];
                ++aK;
            } else if (aI[aK] < bI[bK]) {
                // a is behind b, advance a
                outputX[outputK] = aX[aK];
                outputI[outputK] = aI[aK];
                ++aK;
            } else if (bI[bK] < aI[aK]) {
                // b is behind a, advance b
                outputX[outputK] = bX[bK];
                outputI[outputK] = bI[bK];
                ++bK;
            } else {
                // they are at the same row, advance both
                outputX[outputK] = aX[aK] + bX[bK];
                outputI[outputK] = bI[bK];
                ++aK;
                ++bK;
            }
            ++outputK;
        }
    }

    S4 output;
    if (a.inherits("dsCMatrix")) {
        output = S4("dsCMatrix");
        output.slot("uplo") = a.slot("uplo");
    } else {
        output = S4("dgCMatrix");
    }
    output.slot("i") = outputI_;
    output.slot("x") = outputX_;
    output.slot("p") = outputP_;
    output.slot("Dim") = IntegerVector({ nRows, nColumns });
    output.slot("factors") = List();
    output.slot("Dimnames") = List({ R_NilValue, R_NilValue });

    return output;
}

// Fast kronecker with general sparse matrices
// [[Rcpp::export(name=".fast_kronecker", rng=false)]]
S4 fast_kronecker(S4 a, S4 b) {
    IntegerVector aDims(a.slot("Dim"));
    IntegerVector aI(a.slot("i"));
    IntegerVector aP(a.slot("p"));
    NumericVector aX(a.slot("x"));

    IntegerVector bDims(b.slot("Dim"));
    IntegerVector bI(b.slot("i"));
    IntegerVector bP(b.slot("p"));
    NumericVector bX(b.slot("x"));

    int nRows = aDims[0] * bDims[0];
    int nColumns = aDims[1] * bDims[1];
    NumericVector outputX(aX.size() * bX.size());
    IntegerVector outputI(aX.size() * bX.size());
    IntegerVector outputP(nColumns + 1);

    int outputK = 0;
    outputP[0] = 0;
    for (int outputCol = 0; outputCol < nColumns; ++outputCol) {
        int aCol = outputCol / bDims[1];
        int bCol = outputCol % bDims[1];

        // Each row in this column of a
        for (int aK = aP[aCol]; aK < aP[aCol + 1]; ++aK) {
            int aRow = aI[aK];
            // Each row in this column of b
            for (int bK = bP[bCol]; bK < bP[bCol + 1]; ++bK) {
                int bRow = bI[bK];

                outputX[outputK] = aX[aK] * bX[bK];
                outputI[outputK] = aRow * bDims[0] + bRow;
                ++outputK;
            }
        }

        outputP[outputCol + 1] = outputK;
    }

    S4 output("dgCMatrix");
    output.slot("Dim") = IntegerVector({ nRows, nColumns });
    output.slot("x") = outputX;
    output.slot("i") = outputI;
    output.slot("p") = outputP;
    output.slot("factors") = List();
    output.slot("Dimnames") = List({ R_NilValue, R_NilValue });

    return output;
}

// Count the number of non-zero elements on the diagonal
int nnzero_diag(S4 x) {
    IntegerVector xDims(x.slot("Dim"));
    IntegerVector xI(x.slot("i"));
    IntegerVector xP(x.slot("p"));

    int count = 0;
    for (int xCol = 0; xCol < xDims[1]; ++xCol) {
        if (xP[xCol + 1] == xP[xCol]) continue;
        if (xI[xP[xCol + 1] - 1] == xCol) ++count;
    }
    return count;
}

// Fast kronecker when a and b are symmetrical
// [[Rcpp::export(name=".fast_kronecker_sym", rng=false)]]
S4 fast_kronecker_sym(S4 a, S4 b) {
    std::string aUplo = a.slot("uplo");
    std::string bUplo = b.slot("uplo");
    if (aUplo != "U" || bUplo != "U") {
        Rcpp::stop("Both a and b must use upper triangle representation");
    }

    IntegerVector aDims(a.slot("Dim"));
    IntegerVector aI(a.slot("i"));
    IntegerVector aP(a.slot("p"));
    NumericVector aX(a.slot("x"));

    IntegerVector bDims(b.slot("Dim"));
    IntegerVector bI(b.slot("i"));
    IntegerVector bP(b.slot("p"));
    NumericVector bX(b.slot("x"));

    int nnzA = aX.size();
    int nnzB = bX.size();
    int nnzADiag = nnzero_diag(a);
    int nnzBDiag = nnzero_diag(b);

    // Count the number of non-zero entries in each row of B
    std::vector<int> nnzBelowDiagB(bDims[0]);
    for (int bCol = 0; bCol < bDims[1]; ++bCol) {
        for (int bK = bP[bCol]; bK < bP[bCol + 1]; ++bK) {
            if (bCol == bI[bK]) continue;
            ++nnzBelowDiagB[bI[bK]];
        }
    }

    int nRows = aDims[0] * bDims[0];
    int nColumns = aDims[1] * bDims[1];
    int outputSize = nnzADiag * nnzB + (nnzA - nnzADiag) * (nnzBDiag + 2 * (nnzB - nnzBDiag));

    NumericVector outputX(outputSize);
    IntegerVector outputI(outputSize);
    IntegerVector outputP(nColumns + 1);

    // Compute outputP
    outputP[0] = 0;
    for (int outputCol = 0; outputCol < nColumns; ++outputCol) {
        int aCol = outputCol / bDims[1];
        int nnzACol = aP[aCol + 1] - aP[aCol];

        int n;
        if (nnzACol == 0) {
            n = 0;
        } else {
            int nnzAColDiag;
            int nnzAColAboveDiag;
            int lastARow = aI[aP[aCol + 1] - 1];
            if (lastARow == aCol) {
                nnzAColDiag = 1;
                nnzAColAboveDiag = aP[aCol + 1] - aP[aCol] - 1;
            } else {
                nnzAColDiag = 0;
                nnzAColAboveDiag = aP[aCol + 1] - aP[aCol];
            }

            int bCol = outputCol % bDims[1];
            int nnzBColUpper = bP[bCol + 1] - bP[bCol];
            int nnzBColLower = nnzBelowDiagB[bCol];

            n = nnzAColDiag * nnzBColUpper + nnzAColAboveDiag * (nnzBColUpper + nnzBColLower);
        }

        outputP[outputCol + 1] = outputP[outputCol] + n;
    }

    std::vector<int> bBelowDiagSoFar(bDims[0]);
    for (int aCol = 0; aCol < aDims[1]; ++aCol) {
        for (int aK = aP[aCol]; aK < aP[aCol + 1]; ++aK) {
            int aRow = aI[aK];
            std::fill(bBelowDiagSoFar.begin(), bBelowDiagSoFar.end(), 0);
            for (int bCol = 0; bCol < bDims[1]; ++bCol) {
                for (int bK = bP[bCol]; bK < bP[bCol + 1]; ++bK) {
                    int bRow = bI[bK];

                    double value = aX[aK] * bX[bK];

                    int outputCol = aCol * bDims[1] + bCol;
                    int outputRow = aRow * bDims[0] + bRow;
                    int outputK = (
                        // Start of column
                        outputP[outputCol]
                        // Number of rows before it
                        + (aK - aP[aCol]) * (bP[bCol + 1] - bP[bCol] + nnzBelowDiagB[bCol])
                        // How far into this row of b so far
                        + bK - bP[bCol]
                    );
                    outputX[outputK] = value;
                    outputI[outputK] = outputRow;

                    if (aCol != aRow && bCol != bRow) {
                        // Off a diagonal, so need to reflect b
                        outputCol = aCol * bDims[1] + bRow;
                        outputRow = aRow * bDims[0] + bCol;
                        outputK = (
                            // Start of column
                            outputP[outputCol]
                            // Number of rows before it
                            + (aK - aP[aCol]) * (bP[bRow + 1] - bP[bRow] + nnzBelowDiagB[bRow])
                            // Offset into the lower triangle
                            + bP[bRow + 1] - bP[bRow]
                            // How far into this col of b so far
                            + bBelowDiagSoFar[bRow]
                        );
                        outputX[outputK] = value;
                        outputI[outputK] = outputRow;

                        ++bBelowDiagSoFar[bRow];
                    }
                }

            }
        }
    }

    S4 output("dsCMatrix");
    output.slot("Dim") = IntegerVector({ nRows, nColumns });
    output.slot("x") = outputX;
    output.slot("i") = outputI;
    output.slot("p") = outputP;
    output.slot("uplo") = "U";
    output.slot("factors") = List();
    output.slot("Dimnames") = List({ R_NilValue, R_NilValue });

    return output;
}
