#include <library.h>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>

// build matrix:
void build_problem(std::vector<Eigen::Triplet<double>> &tripletList, const std::vector<scalar>& r, int N, double dr)
{

    tripletList.reserve(3 * N);

    // loop over interior nodes and fill the matrix:
    for (int i = 1; i < N - 1; ++i)
    {
        tripletList.push_back({i, i, -2.0});
        tripletList.push_back({i, i - 1, 1.0 - (0.5*dr)/r[i]});
        tripletList.push_back({i, i + 1, 1.0 + (0.5*dr)/r[i]});
    }

    // compute coefficients for left boundary nodes:
    tripletList.push_back({0, 0, -1.5*dr});
    tripletList.push_back({0, 1, 2*dr});
    tripletList.push_back({0, 2, -0.5*dr});

    // compute coefficients for right boundary nodes:
    tripletList.push_back({N - 1, N - 1, 1.0});


    tripletList.shrink_to_fit();

}


int main()
{
    // read simulation parameters:
    auto parameters = read_input("input");
    scalar S = parameters[0];
    int N = static_cast<int>(parameters[1]);

    // compute grid size:
    scalar dr = 1.0 / static_cast<double>(N - 1);

    // create grid points:
    std::vector<scalar> r = linspace(0.0, 1.0, N);

    // build A matrix:
    Eigen::SparseMatrix<double, Eigen::RowMajor> A(N, N);
    std::vector<Eigen::Triplet<double>> tripletList;
    build_problem(tripletList, r, N, dr);
    A.setFromTriplets(tripletList.begin(), tripletList.end());

    // build B vector
    Eigen::VectorXd B = Eigen::VectorXd::Constant(N, -S*dr*dr);

    // Boundary condition:
    B(0) = 0.0;
    B(N-1)= 1.0;

    // solve AX = B:
    Eigen::SparseLU<Eigen::SparseMatrix<scalar, Eigen::RowMajor>> solver;
    solver.analyzePattern(A);
    solver.factorize(A);
    Eigen::VectorXd T(N);
    T = solver.solve(B);
    assert(solver.info() == Eigen::Success);


    // write x vs T for x:
    write_data write(std::to_string(S) + ".dat");

    for (int i = 0; i < N; ++i)
        write(r[i], T[i]);


    return 0;
}
