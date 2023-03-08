// write data to be plotted in python or paraview:
#pragma once

#include <armadillo>
#include <fstream>


// write data to be plotted using matplotlib library:
void write_matplotlib(const arma::vec &x, const arma::vec &y, const arma::vec &sol)
{
    int nx = x.size(), ny = y.size();
    auto map = [nx](int i, int j)
    { return nx * j + i; };

    std::ofstream file("sol.dat");
    file << "x,y,sol\n";
    file.flags(std::ios::dec | std::ios::scientific);
    file.precision(16);

    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            file << x(i) << " " << y(j) << " " << sol(map(i, j)) << "\n";

    file.close();
}

// write data in vtk format:
void write_vtk(const arma::vec &x, const arma::vec &y, const arma::vec &sol)
{
    int nx = x.size(), ny = y.size();
    auto map = [nx](int i, int j)
    { return nx * j + i; };

    std::ofstream file("sol.vtk");
    file.flags(std::ios::dec | std::ios::scientific);
    file.precision(16);

    file << "# vtk DataFile Version 3.0\n";
    file << "Structured Grid Example\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_GRID\n";
    file << "DIMENSIONS " << nx << " " << ny << " 1\n";
    file << "POINTS " << nx * ny << " float\n";

    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            file << x(i) << " " << y(j) << " 0\n";

    file << "POINT_DATA " << nx * ny << "\n";
    file << "SCALARS Solution float 1\n";
    file << "LOOKUP_TABLE default\n";

    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            file << sol(map(i, j)) << "\n";

    file.close();
}