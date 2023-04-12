#pragma once

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

    // Open file for writing
    std::ofstream file(filename.str(), std::ios::out);
    file.flags(std::ios::dec | std::ios::scientific);
    file.precision(16);

    if (!file)
        throw std::runtime_error("Cannot open output file!");

    // Write numerical and exact solutions to file
    int n = x.size();
    for (int i = 0; i < n; i++)
        file << time << " " << x[i] << " " << u[i] << " " << uexact[i] << std::endl;

    // Close file
    file.close();
}

// compute average error:
double average_error(const std::vector<double>& u, const std::vector<double>& uexact)
{
    int N = u.size();
    double error = 0.0;

    for (int i = 0; i < N; ++i)
        error += fabs(u[i] - uexact[i]);

    return error/N;
}
