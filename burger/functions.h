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
inline std::vector<double> linspace_periodic(double xmin, double xmax, int n)
{
    double h = (xmax - xmin) / static_cast<double>(n);
    std::vector<double> x(n);

    for (int i = 0; i < n; ++i)
        x[i] = xmin + i * h;

    return x;
}

void write_matplotlib(int iteration, double time, const std::vector<double> &x, const std::vector<double> &u)
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
        file << time << " " << x[i] << " " << u[i] << std::endl;

    // Close file
    file.close();
}
