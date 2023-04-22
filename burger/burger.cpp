// solve 1d burger's equation on a periodic domain using the FTCS scheme:

#include "functions.h"

// burger update using upwinding:
void burger_upwind(const std::vector<double>& u, std::vector<double>& un, double rd, int n, double h, double dt)
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

}

// burger update using upwinding:
void burger_central(const std::vector<double>& u, std::vector<double>& un, double rd, int n, double h, double dt)
{
    // update solution:
        for (int i = 1; i < n - 1; ++i)
            un[i] = u[i] + rd * (u[i + 1] - 2.0 * u[i] + u[i - 1])
                         - (0.5*dt/h)*(u[i + 1] - u[i - 1]);

        // periodic boundary condition:
        un[0] = u[0] + rd * (u[1] - 2.0 * u[0] + u[n - 1])
                     - (0.5*dt/h)*(u[1] - u[n - 1]);



        un[n - 1] = u[n - 1] + rd * (u[0] - 2.0 * u[n - 1] + u[n - 2])
                             - (0.5*dt/h)*(u[0] - u[n - 2]);
                           
}



int main()
{
    // discretize x:
    const double xmin = 0.0, xmax = 1.0;
    const int n = 1024;
    const double h = (xmax - xmin) / static_cast<double>(n);
    const std::vector<double> x = linspace_periodic(xmin, xmax, n);

    // discretize t:
    double t = 0.0;
    const double tf = 1.0;
    const double dt = 0.0004;


    // stability parameter
    const double alpha = 0.000;
    const double rd = alpha*dt/(h*h);


    // exact solution:
    auto ic = [](double x)
    { return sin(2.0 * M_PI * x); };
    //{ return sin(4.0 * M_PI * x) + sin(6.0 * M_PI * x) + sin(10.0 * M_PI * x); };
   
    

    // initial condition:
    std::vector<double> u(n), un(n);
    std::transform(x.begin(), x.end(), u.begin(), ic);
    write_matplotlib(0, t, x, u);
    
    // evolve in time:
    int it = 0;
    while (t < tf)
    {
        burger_upwind(u, un, rd, n, h, dt);

        // update time:
        t += dt;
        it++;

        // swap:
        u.swap(un);

        // write data:
        if (it % 100 == 0) write_matplotlib(it, t, x, u);

    }

    // compute error:   
    std::ofstream file("error.dat", std::ios::app);
    file.flags(std::ios::dec | std::ios::scientific);
    file.precision(16);
    

    return 0;
}
