// solve diffusion equation using BTCS on a periodic domain. Use Jacobi method to solve linear system.
#include "functions.h"

int jacobi(std::vector<double> &un, const std::vector<double> &u, double rd)
{
    int n = un.size();
    int iter = 0;
    double coeff = (1.0 + 2.0 * rd);
    double invcoeff = 1.0 / coeff;
    std::vector<double> uiter(u);
    std::vector<double> residue(n);

    // jacobi iteration:
    do
    {
        // update initerior nodal values:
        for (int i = 1; i < n - 1; ++i)
        {
            un[i] = (rd * uiter[i + 1] + rd * uiter[i - 1] + u[i]) * invcoeff;
            residue[i] = fabs(u[i] - coeff * un[i] + rd * un[i + 1] + rd * un[i - 1]);
        }

        // periodic boundary:

        // update left value:
        un[0] = (rd * uiter[1] + rd * uiter[n - 1] + u[0]) * invcoeff;
        residue[0] = fabs(u[0] - coeff * un[0] + rd * un[1] + rd * un[n - 1]);

        // update right value:
        un[n - 1] = (rd * uiter[0] + rd * uiter[n - 2] + u[n - 1]) * invcoeff;
        residue[n - 1] = fabs(u[n - 1] - coeff * un[n - 1] + rd * un[0] + rd * un[n - 2]);

        uiter.swap(un);

        iter++;

    } while (*std::max_element(residue.begin(), residue.end()) > 1.0e-9);

    return iter;
}

int main()
{
    // discretize x:
    const double xmin = 0.0, xmax = 2.0 * M_PI;
    const int n = 256;
    const double h = (xmax - xmin) / static_cast<double>(n);
    const std::vector<double> x = linspace_periodic(xmin, xmax, n);

    // diffusion coefficient:
    const double alpha = 1.0;

    // stability parameter
    const double rd = 0.5;

    // discretize t:
    double t = 0.0;
    const double tf = 0.4;
    const double dt = (rd * h * h) / alpha;

    // exact solution:
    auto f = [&t, alpha](double x)
    { return exp(-4.0 * alpha * t) * sin(2.0 * x); };

    // initial condition:
    std::vector<double> u(n), un(n), uexact(n);
    std::transform(x.begin(), x.end(), u.begin(), f);
    std::transform(x.begin(), x.end(), uexact.begin(), f);
    write_matplotlib(0, 0, x, u, uexact);

    // evolve in time:
    int it = 0;
    while (t < tf)
    {
        // solve the linear system using the jacobi method:
        jacobi(un, u, rd);

        // update time:
        t += dt;
        it++;

        // swap:
        u.swap(un);

        // exact solution:
        std::transform(x.begin(), x.end(), uexact.begin(), f);

        // write data:
        if (it % 20 == 0)
            write_matplotlib(it, t, x, u, uexact);
    }

    // compute error:
    std::ofstream file("error.dat", std::ios::app);
    file.flags(std::ios::dec | std::ios::scientific);
    file.precision(16);

    file << n << "\t" << average_error(u, uexact) << "\n";

    return 0;
}