// solve 1d diffusion equation on a periodic domain using the FTCS scheme:

#include "functions.h"


int main()
{
    // discretize x:
    const double xmin = 0.0, xmax = 2.0 * M_PI;
    const int n = 256;
    const double h = (xmax - xmin) / static_cast<double>(n - 1);
    const std::vector<double> x = linspace(xmin, xmax, n);

    // diffusion coefficient:
    const double alpha = 1.0;

    // stability parameter
    const double rd = 0.166667;

    // discretize t:
    double t = 0.0;
    const double tf = 1;
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
        // update solution:
        for (int i = 1; i < n - 1; ++i)
            un[i] = u[i] + rd * (u[i + 1] - 2.0 * u[i] + u[i - 1]);

        // TODO introduce periodic boundary condition:
        //! code still works because of zero values at the boundary!

        // update time:
        t += dt;
        it++;

        // swap:
        u.swap(un);

        // exact solution:
        std::transform(x.begin(), x.end(), uexact.begin(), f);

        // write data:
        //if (it % 20 == 0) write_matplotlib(it, t, x, u, uexact);

    }

    // compute error:   
    std::ofstream file("error.dat", std::ios::app);
    file.flags(std::ios::dec | std::ios::scientific);
    file.precision(16);
    
    file << n << "\t" << average_error(u, uexact) << "\n";
    

    return 0;
}
