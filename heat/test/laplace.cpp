// solve poisson equation for Dirchlet boundary conditions:

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <unsupported/Eigen/SparseExtra>

#include <iostream>
#include <fstream>
#include <cassert>

#include "plot.h"

using namespace Eigen;
using namespace std;

using SpMat = SparseMatrix<double, RowMajor>;
using Vec = VectorXd;

void build_matrix(SpMat &A, const int nx, const int ny)
{
    // grid index to vec index map:
    auto map = [nx](int i, int j)
    { return nx * j + i; };

    // create tripletList:
    vector<Triplet<double>> tripleList;
    tripleList.reserve(nx * ny * 5);

    // interior nodes:
    for (int i = 1; i < nx - 1; ++i)
    {
        for (int j = 1; j < ny - 1; ++j)
        {
            tripleList.emplace_back(map(i, j), map(i, j), -4.0);
            tripleList.emplace_back(map(i, j), map(i - 1, j), 1.0);
            tripleList.emplace_back(map(i, j), map(i + 1, j), 1.0);
            tripleList.emplace_back(map(i, j), map(i, j - 1), 1.0);
            tripleList.emplace_back(map(i, j), map(i, j + 1), 1.0);
        }
    }

    // for boundary nodes set diagonal elements = 1:
    // top and bottom boundary: dirichlet boundary
    for (int i = 1; i < nx - 1; ++i)
    {
        int j = 0;
        tripleList.emplace_back(map(i, j), map(i, j), 1.0);
        j = ny - 1;
        tripleList.emplace_back(map(i, j), map(i, j), 1.0);
    }

    // left and right boundary: dirichlet boundary
    for (int j = 0; j < ny; ++j)
    {
        int i = 0;
        tripleList.emplace_back(map(i, j), map(i, j), 1.0);
        i = nx - 1;
        tripleList.emplace_back(map(i, j), map(i, j), 1.0);
    }

    A.setFromTriplets(tripleList.begin(), tripleList.end());
}

int main()
{
    // exact solution:
    auto uexact = [](const double &x, const double &y)
    { return 1.0; };

    // create x and y grid pts:
    const int nx = 3;
    const int ny = 3;
    const double xmin = 0.0, xmax = 1.0;
    const double ymin = 0.0, ymax = 1.0;

    auto x = ArrayXd::LinSpaced(nx, xmin, xmax);
    auto y = ArrayXd::LinSpaced(ny, ymin, ymax);

    // compute grid size:
    const double h = (xmax - xmin) / static_cast<double>(nx - 1);

    // grid index to vec index map:
    auto map = [nx](int i, int j)
    { return nx * j + i; };

    // size of solution vector:
    const int N = nx * ny;

    // initialize A matirx and B vector:
    SpMat A(N, N);
    Vec B = Vec::Zero(N);

    // build A matrix:
    build_matrix(A, nx, ny);

    std::cout << A << std::endl;

    // build vector B:

    // top and bottom:
    for (int i = 1; i < nx - 1; ++i)
    {
        int j = 0;
        B(map(i, j)) = uexact(x(i), y(j));
        j = ny - 1;
        B(map(i, j)) = uexact(x(i), y(j));
    }

    // left and right:
    for (int j = 0; j < ny; ++j)
    {
        int i = 0;
        B(map(i, j)) = uexact(x(i), y(j));
        i = nx - 1;
        B(map(i, j)) = uexact(x(i), y(j));
    }

    std::cout << B << "\n\n";

    // solve AX = B:
    SparseLU<decltype(A)> solver;
    solver.compute(A);
    Vec X = solver.solve(B);


    std::cout << X << std::endl;

    // write L2 and Linfty error to a file:
    Vec Xexact(N);
    for (int i = 0; i < nx; ++i)
        for (int j = 0; j < ny; ++j)
            Xexact(map(i, j)) = uexact(x(i), y(j));

    Vec Error = Xexact - X;

    std::ofstream file("error.dat", std::ios::app);
    file.flags(std::ios::dec | std::ios::scientific);

    file << h << "\t" << h * Error.norm() << "\t" << Error.cwiseAbs().maxCoeff() << std::endl;
    file.close();

    // postprocessing:
    write_matplotlib(x, y, X);

    return 0;
}
