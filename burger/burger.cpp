// solve 1d burger's equation on a periodic domain using the FTCS scheme:

#include "functions.h"

int main()
{
    // discretize x:
    const double xmin = 0.0, xmax = 1.0;
    const int n = 64;
    const double h = (xmax - xmin) / static_cast<double>(n);
    const std::vector<double> x = linspace_periodic(xmin, xmax, n);

    // stability parameter
    const double rd = 0.166667;

    // discretize t:
    double t = 0.0;
    const double tf = 0.075;
    const double dt = 0.0004;

    // exact solution:
    auto ic = [](double x)
    { return sin(4.0 * M_PI * x) + sin(6.0 * M_PI * x) + sin(10.0 * M_PI * x); };

    // initial condition:
    std::vector<double> u(n), un(n);
    std::transform(x.begin(), x.end(), u.begin(), ic);
    
    // evolve in time:
    int it = 0;
    while (t < tf)
    {
        // update solution:
        for (int i = 1; i < n - 1; ++i)
            un[i] = u[i] + rd * (u[i + 1] - 2.0 * u[i] + u[i - 1])
                         - (dt/h)*std::max(u[i], 0.0)*(u[i] - u[i - 1])
                         - (dt/h)*std::min(u[i], 0.0)*(u[i + 1] - u[i]);

        // periodic boundary condition:
        un[0] = u[0] + rd * (u[1] - 2.0 * u[0] + u[n - 1])
                     - (dt/h)*std::max(u[0], 0.0)*(u[0] - u[n - 1])
                     - (dt/h)*std::min(u[0], 0.0)*(u[1] - u[0]);



        un[n - 1] = u[n - 1] + rd * (u[0] - 2.0 * u[n - 1] + u[n - 2])
                             - (dt/h)*std::max(u[n - 1], 0.0)*(u[n - 1] - u[n - 2])
                             - (dt/h)*std::min(u[n - 1], 0.0)*(u[0] - u[n - 1]);


        // update time:
        t += dt;
        it++;

        // swap:
        u.swap(un);

        // write data:
        if (it % 10 == 0) write_matplotlib(it, t, x, u);

    }

    // compute error:   
    std::ofstream file("error.dat", std::ios::app);
    file.flags(std::ios::dec | std::ios::scientific);
    file.precision(16);
    

    return 0;
}
