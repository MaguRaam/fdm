// solve 1d diffusion equation on a periodic domain using the FTCS scheme:

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <iomanip>

// create grid points:
inline std::vector<double> linspace(double xmin, double xmax, int n)
{
    double h = (xmax - xmin) / static_cast<double>(n - 1);
    std::vector<double> x(n);

    for (int i = 0; i < n; ++i)
        x[i] = xmin + i * h;

    return x;
}

// write data in vtk format:
void write_vtk(int it, const std::vector<double> &x, const std::vector<double> &u, const std::vector<double> &uexact)
{
    const int n = x.size();
    const std::string filename = "output_" + std::to_string(it) + ".vtk";
    std::ofstream outfile(filename);
    
    if (!outfile) {
    	throw std::runtime_error("Cannot open output file " + filename + "!");
		}

    // write vtk header:
    outfile << "# vtk DataFile Version 3.0\n";
    outfile << "vtk output\n";
    outfile << "ASCII\n";
    outfile << "DATASET POLYDATA\n";

    // write grid points:
    outfile << "POINTS " << n << " float\n";
    for (int i = 0; i < n; ++i)
        outfile << x[i] << " 0.0 0.0\n";

    // write solution vectors:
    outfile << "POINT_DATA " << n << "\n";
    outfile << "SCALARS u float\n";
    outfile << "LOOKUP_TABLE default\n";
    for (int i = 0; i < n; ++i)
        outfile << u[i] << "\n";
    outfile << "SCALARS uexact float\n";
    outfile << "LOOKUP_TABLE default\n";
    for (int i = 0; i < n; ++i)
        outfile << uexact[i] << "\n";
}

void write_matplotlib(int iteration, double time, const std::vector<double> &x, const std::vector<double> &u,
                      const std::vector<double> &uexact)
{
    // Create filename based on iteration number and time
    std::ostringstream filename;
    filename << "data/data_" << std::setw(5) << std::setfill('0') << iteration << ".dat";

		if (!outfile) {
    	throw std::runtime_error("Cannot open output file " + filename + "!");
		}

    // Open file for writing
    std::ofstream file(filename.str(), std::ios::out);
    file.flags(std::ios::dec | std::ios::scientific);
    file.precision(16);

    // Write numerical and exact solutions to file
    int n = x.size();
    for (int i = 0; i < n; i++)
        file << time << " " << x[i] << " " << u[i] << " " << uexact[i] << std::endl;

    // Close file
    file.close();
}

int main()
{
    // discretize x:
    const double xmin = 0.0, xmax = 2.0 * M_PI;
    const int n = 64;
    const double h = (xmax - xmin) / static_cast<double>(n - 1);
    const std::vector<double> x = linspace(xmin, xmax, n);

    // diffusion coefficient:
    const double alpha = 1.0;

    // stability parameter
    const double rd = 0.5;

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
        if (it % 20 == 0)
            write_matplotlib(it, t, x, u, uexact);
    }

    return 0;
}
