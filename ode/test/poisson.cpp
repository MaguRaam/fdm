#include <iostream>
#include <armadillo>


int main()
{
    // create grid points
    const int N = 64000;
    const arma::vec x = arma::linspace(0.0, 1.0, N);

    // compute grid size:
    const double h = 1.0 / static_cast<double>(N - 1);

    // create exact solution:
    auto exact = [](double x)
    { return sin(M_PI * x); };

    // create rhs function:
    auto rhs = [h](double x)
    { return -M_PI * M_PI * sin(M_PI * x) * h * h; };

    // create exact solution vector:
    arma::vec u_exact(N);
    std::transform(x.begin(), x.end(), u_exact.begin(), exact);

    // create vector b:
    arma::vec b(N);
    std::transform(x.begin(), x.end(), b.begin(), rhs);

    // boundary consditions;
    b(0) = b(N - 1) = 0.0;

    // create sparse matrix A:
    arma::sp_mat A(N, N);

    for (int i = 1; i < N - 1; ++i)
    {
        A(i, i) = -2;
        A(i, i - 1) = +1;
        A(i, i + 1) = +1;
    }

    // boundary conditions:
    A(0, 0) = A(N - 1, N - 1) = 1.0;

    //solve:
    arma::vec u = arma::spsolve(A,b, "superlu");

    // write data to file:
    std::fstream file("sol.dat", std::ios::out);
    file.flags(std::ios::dec | std::ios::scientific);
    file.precision(16);

    for (int i = 0; i < N; ++i)
        file << x(i) << "\t" << u(i) << "\t" << u_exact(i) << "\n";
    file.close();

    // compute error:
    file.open("error.dat", std::ios::app);
    file << h << "\t" << arma::abs(u - u_exact).max() << "\n";


    return 0;
}
