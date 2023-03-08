// Solve Laplace equation using Direct Solver 

#include <iostream>
#include <armadillo>
#include "plot.h"

void build_matrix(arma::sp_mat &A, const int nx, const int ny, const double h)
{
    // grid index to vec index map:
    auto map = [nx](int i, int j)
    { return nx * j + i; };

    // interior nodes:
    for (int i = 1; i < nx - 1; ++i)
    {
        for (int j = 1; j < ny - 1; ++j)
        {
            A(map(i, j), map(i, j)) = -4.0;
            A(map(i, j), map(i - 1, j)) = 1.0;
            A(map(i, j), map(i + 1, j)) = 1.0;
            A(map(i, j), map(i, j + 1)) = 1.0;
            A(map(i, j), map(i, j - 1)) = 1.0;
        }
    }

    // boundary nodes:

    // left and right boundary is dirichlet, set diagonal element = 1
    for (int j = 0; j < ny; ++j)
    {
        int i = 0;
        A(map(i, j), map(i, j)) = 1.0;
        i = nx - 1;
        A(map(i, j), map(i, j)) = 1.0;
    }


    // top and bottom boundary is neumann, solve laplace equation using ghost point:
    for (int i = 1; i < nx - 1; ++i)
    {
        // bottom:
        int j = 0;
        A(map(i, j), map(i, j)) = -4.0;
        A(map(i, j), map(i - 1, j)) = 1.0;
        A(map(i, j), map(i + 1, j)) = 1.0;
        A(map(i, j), map(i, j + 1)) = 2.0;

        // top:
        j = ny - 1;
        A(map(i, j), map(i, j)) = 2*h -4.0;
        A(map(i, j), map(i - 1, j)) = 1.0;
        A(map(i, j), map(i + 1, j)) = 1.0;
        A(map(i, j), map(i, j - 1)) = 2.0;
    }

}

int main()
{
    // create x and y grid pts:
    const int nx = 128;
    const int ny = 256;
    const double xmin = 0.0, xmax = 1.0;
    const double ymin = 0.0, ymax = 2.0;
    const arma::vec x = arma::linspace(xmin, xmax, nx);
    const arma::vec y = arma::linspace(ymin, ymax, ny);

    // compute grid size:
    const double h = (xmax - xmin) / static_cast<double>(nx - 1);

    // size of solution vector:
    const int N = nx * ny;

    // grid index to vec index map:
    auto map = [nx](int i, int j)
    { return nx * j + i; };

    // create sparse matrix A:
    arma::sp_mat A(N, N);
    arma::vec B(N);

    // build matrix:
    build_matrix(A, nx, ny, h);

    // build vector B:

    // left and right:
    for (int j = 0; j < ny; ++j)
    {
        int i = 0;
        B(map(i, j)) = 30.0;
        i = nx - 1;
        B(map(i, j)) = 60.0;
    }

    // top:
    for (int i = 1; i < nx - 1; ++i)
    {
        int j = ny - 1;
        B(map(i, j)) = 120*h;
    }

    // solve AX = B:
    arma::vec X = arma::spsolve(A, B, "superlu");

    // write data:
    write_matplotlib(x, y, X);

    return 0;
}
