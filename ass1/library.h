#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cassert>
#include <sstream>
#include <cmath>
#include <numeric>

// ODE solvers:

//!Note: The ODE solvers take input as a function of the form f(y, t):

// update state using explicit Euler method:
template <typename F, typename State, typename Scalar>
inline void euler_update(F f, State &y, const Scalar t, const Scalar dt)
{
    y = y + dt*f(y, t);
}


// update state using RK2 method:
template <typename F, typename State, typename Scalar>
inline void rk2_update(F f, State &y, const Scalar t, const Scalar dt)
{

    State k1 = f(y, t);
    State k2 = f(y + 0.5 * dt* k1, t + 0.5*dt);
    
    y = y + dt*k2;
}

// update state using RK4 method:
template <typename F, typename State, typename Scalar>
inline void rk4_update(F f, State &y, const Scalar t, const Scalar dt)
{

    State k1 = f(y, t);
    State k2 = f(y + 0.5 * dt * k1, t + 0.5 * dt);
    State k3 = f(y + 0.5 * dt * k2, t + 0.5 * dt);
    State k4 = f(y + dt * k3, t + dt);

    y = y + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
}

// update state using 2nd order Adams-Bashforth scheme:
template <typename F, typename State, typename Scalar>
inline State ab2_update(F f, const State &ycurrent, const State & yold, const Scalar tcurrent, const Scalar dt)
{
    Scalar told = tcurrent - dt;
    return ycurrent + 0.5*dt*( 3.0*f(ycurrent, tcurrent) - f(yold, told));
}


// state type:
template <typename T>
struct state
{
    T x, v;
};

// add two states: s1 + s2
template <typename T>
state<T> operator+(const state<T> &s1, const state<T> &s2)
{
    return {s1.x + s2.x, s1.v + s2.v};
}

// scalar multiplication with state: a*s
template <typename T, typename W>
state<T> operator*(const W &a, const state<T> &s)
{
    return {a * s.x, a * s.v};
}

// scalar multiplication with state: s*a
template <typename T, typename W>
state<T> operator*(const state<T> &s, const W &a)
{
    return {a * s.x, a * s.v};
}

// scalar type:
using scalar = double;

// IO:

// read simulation parameters from input file:
std::vector<scalar> read_input(std::string input_filename)
{
    // open parameter file:
    std::ifstream infile(input_filename);
    assert(infile.is_open() && "input file doesn't exist");

    // vector to store all the input parameters:
    std::vector<scalar> parameters;

    // read file line by line:
    std::string line;
    scalar value;
    while (std::getline(infile, line))
    {
        // skip comments and empty lines:
        if (line.find("!") != std::string::npos || line == "")
            continue;

        // push back values to parameter vector:
        std::istringstream iss(line);
        iss >> value;
        parameters.push_back(value);
    }

    return parameters;
}

// functor to write arbirtary no of columns to a file:
class write_data
{
public:
    write_data(std::string filename) : file(filename)
    {
        file.flags(std::ios::dec | std::ios::scientific);
        file.precision(16);
    }

    void operator()(){file << "\n";}
    
    // variadic function call operator to write arbirtary no of columns to file:
    template <typename T, typename... Types>
    void operator()(T var1, Types... var2)
    {
        file << var1 << "\t";
        (*this)(var2...);
    }

private:
    std::ofstream file;
};

// create a grid of points (similar to linspace function in python)
inline std::vector<scalar> linspace(scalar start, scalar stop, int n)
{
    std::vector<scalar> result(n);
    scalar step = (stop - start) / static_cast<scalar>(n - 1);
    for (int i = 0; i < n; ++i)
        result[i] = start + i * step;
    return result;
}



