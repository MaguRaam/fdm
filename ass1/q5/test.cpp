#include <library.h>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>

int main()
{

    // no of grid points:
    int N = 3200;

    // compute grid size:
    scalar h = 1.0 / static_cast<scalar>(N - 1);

    // create grid points:
    std::vector<scalar> x = linspace(0.0, 1.0, N);

    // exact solution:
    auto uexact = [](scalar xi)
    { return sin(M_PI * xi); };

    // rhs function:
    auto rhs = [h](scalar xi)
    { return -h * h * M_PI * M_PI * sin(M_PI * xi); };

    // create A, X, Xexact and B:
    Eigen::SparseMatrix<scalar, Eigen::RowMajor> A(N, N);
    Eigen::VectorXd X = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd Xexact = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd B = Eigen::VectorXd::Zero(N);

    // build A matrix:
    std::vector<Eigen::Triplet<scalar>> tripletList;
    tripletList.reserve(3 * N);

    // loop over interior nodes and fill the matrix:
    for (int i = 1; i < N - 1; ++i)
    {
        tripletList.push_back({i, i, -2.0});
        tripletList.push_back({i, i - 1, 1.0});
        tripletList.push_back({i, i + 1, 1.0});
    }

    // boundary nodes:
    tripletList.push_back({0, 0, 1.0});
    tripletList.push_back({N - 1, N - 1, 1.0});

    // fill sparse matrix:
    tripletList.shrink_to_fit();
    A.setFromTriplets(tripletList.begin(), tripletList.end());

    // build B vector:
    std::transform(x.begin(), x.end(), B.begin(), rhs);

    // Boundary condition:
    B(0) = 0.0;
    B(N - 1) = 0.0;

    // compute Xexact:
    std::transform(x.begin(), x.end(), Xexact.begin(), uexact);

    // solve AX = B:
    Eigen::SparseLU<Eigen::SparseMatrix<scalar, Eigen::RowMajor>> solver;
    solver.analyzePattern(A);
    solver.factorize(A);
    X = solver.solve(B);
    assert(solver.info() == Eigen::Success);

    write_data write("sol.dat");
    for (int i = 0; i < N; ++i)
        write(x[i], X[i], Xexact[i]);

    scalar error = (X - Xexact).cwiseAbs().maxCoeff();
    std::ofstream file("error.dat", std::ios::app);
    file << h << "\t" << error << "\n";
}